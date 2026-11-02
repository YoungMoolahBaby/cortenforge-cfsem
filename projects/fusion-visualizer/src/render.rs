//! A small z-buffer rasterizer: flat-shaded triangle meshes to a 24-bit BMP, for looking at the
//! exported geometry without a 3D viewer.

use cortenforge::mesh::types::IndexedMesh;

pub struct Scene {
    pub items: Vec<(IndexedMesh, [u8; 3])>,
}

/// Render with an orthographic camera looking along `view` (from the eye towards the scene), up
/// roughly +Z, fitting everything into `w x h` pixels.
pub fn render_bmp(scene: &Scene, view: [f64; 3], w: usize, h: usize, path: &std::path::Path) {
    let norm = |v: [f64; 3]| {
        let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        [v[0] / l, v[1] / l, v[2] / l]
    };
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let fwd = norm(view);
    let right = norm(cross(fwd, [0.0, 0.0, 1.0]));
    let up = cross(right, fwd);
    let light = norm([-0.4, -0.6, 0.7]);
    // Screen coordinates (u right, v up, d depth along fwd).
    let proj = |p: [f64; 3]| (dot(p, right), dot(p, up), dot(p, fwd));
    let (mut umin, mut umax, mut vmin, mut vmax) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for (m, _) in &scene.items {
        for p in &m.vertices {
            let (u, v, _) = proj([p.x, p.y, p.z]);
            umin = umin.min(u);
            umax = umax.max(u);
            vmin = vmin.min(v);
            vmax = vmax.max(v);
        }
    }
    let scale = 0.95 * (w as f64 / (umax - umin)).min(h as f64 / (vmax - vmin));
    let (uc, vc) = (0.5 * (umin + umax), 0.5 * (vmin + vmax));
    let to_px = |u: f64, v: f64| {
        (
            0.5 * w as f64 + (u - uc) * scale,
            0.5 * h as f64 - (v - vc) * scale,
        )
    };
    let mut depth = vec![f64::INFINITY; w * h];
    let mut img = vec![[255u8, 255, 255]; w * h];
    for (m, color) in &scene.items {
        for f in &m.faces {
            let p = f.map(|i| {
                let v = m.vertices[i as usize];
                [v.x, v.y, v.z]
            });
            let n = norm(cross(
                [p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]],
                [p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]],
            ));
            if !n[0].is_finite() {
                continue;
            }
            // Two-sided Lambert shading.
            let shade = 0.25 + 0.75 * dot(n, light).abs();
            let c = color.map(|x| (x as f64 * shade).min(255.0) as u8);
            let s: Vec<(f64, f64, f64)> = p
                .iter()
                .map(|&q| {
                    let (u, v, d) = proj(q);
                    let (x, y) = to_px(u, v);
                    (x, y, d)
                })
                .collect();
            let (x0, x1) = (
                s.iter().map(|q| q.0).fold(f64::MAX, f64::min),
                s.iter().map(|q| q.0).fold(f64::MIN, f64::max),
            );
            let (y0, y1) = (
                s.iter().map(|q| q.1).fold(f64::MAX, f64::min),
                s.iter().map(|q| q.1).fold(f64::MIN, f64::max),
            );
            let area =
                (s[1].0 - s[0].0) * (s[2].1 - s[0].1) - (s[2].0 - s[0].0) * (s[1].1 - s[0].1);
            if area.abs() < 1e-12 {
                continue;
            }
            for py in (y0.floor().max(0.0) as usize)..=(y1.ceil().min(h as f64 - 1.0) as usize) {
                for px in (x0.floor().max(0.0) as usize)..=(x1.ceil().min(w as f64 - 1.0) as usize)
                {
                    let (fx, fy) = (px as f64 + 0.5, py as f64 + 0.5);
                    let e = |a: (f64, f64, f64), b: (f64, f64, f64)| {
                        (b.0 - a.0) * (fy - a.1) - (fx - a.0) * (b.1 - a.1)
                    };
                    let w0 = e(s[1], s[2]) / area;
                    let w1 = e(s[2], s[0]) / area;
                    let w2 = 1.0 - w0 - w1;
                    if w0 < 0.0 || w1 < 0.0 || w2 < 0.0 {
                        continue;
                    }
                    let d = w0 * s[0].2 + w1 * s[1].2 + w2 * s[2].2;
                    let k = py * w + px;
                    if d < depth[k] {
                        depth[k] = d;
                        img[k] = c;
                    }
                }
            }
        }
    }
    write_bmp(&img, w, h, path);
}

fn write_bmp(img: &[[u8; 3]], w: usize, h: usize, path: &std::path::Path) {
    let row = (3 * w).div_ceil(4) * 4;
    let size = 54 + row * h;
    let mut b = Vec::with_capacity(size);
    b.extend_from_slice(b"BM");
    b.extend_from_slice(&(size as u32).to_le_bytes());
    b.extend_from_slice(&0u32.to_le_bytes());
    b.extend_from_slice(&54u32.to_le_bytes());
    b.extend_from_slice(&40u32.to_le_bytes());
    b.extend_from_slice(&(w as i32).to_le_bytes());
    b.extend_from_slice(&(h as i32).to_le_bytes());
    b.extend_from_slice(&1u16.to_le_bytes());
    b.extend_from_slice(&24u16.to_le_bytes());
    b.extend_from_slice(&[0u8; 24]);
    for y in (0..h).rev() {
        for x in 0..w {
            let [r, g, bl] = img[y * w + x];
            b.extend_from_slice(&[bl, g, r]);
        }
        b.resize(b.len() + row - 3 * w, 0);
    }
    std::fs::write(path, b).expect("write BMP");
}
