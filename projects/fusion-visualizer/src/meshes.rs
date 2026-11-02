//! `mesh`: the tokamak as CortenForge meshes. Coils are cf-design solids (marching cubes, welded);
//! the vessel wall and flux surfaces are revolved equilibrium contours built as IndexedMesh tori.
//! Each part is checked (watertight, manifold, volume against an exact formula) and written as
//! STL in millimetres; a cut-away render goes to a BMP.

use std::f64::consts::PI;
use std::time::Instant;

use cortenforge::cf_design::Solid;
use cortenforge::mesh::io::save_mesh;
use cortenforge::mesh::repair::{remove_unreferenced_vertices, validate_mesh, weld_vertices};
use cortenforge::mesh::types::{IndexedMesh, Point3};
use nalgebra::{UnitQuaternion, Vector3};

use crate::contour::{contours, inside, perimeter};
use crate::equilibrium::{Coil, Options, solve, summarize};
use crate::gs::Solver;
use crate::render::{Scene, render_bmp};
use crate::shape;
use crate::tf::TfSet;

const MM: f64 = 1000.0;

/// Resample a closed polyline (first point repeated last) to `m` points evenly spaced in arc
/// length, without the repeat.
fn resample(c: &[(f64, f64)], m: usize) -> Vec<(f64, f64)> {
    let total = perimeter(c);
    let mut out = Vec::with_capacity(m);
    let (mut seg, mut acc) = (0, 0.0);
    for k in 0..m {
        let target = total * k as f64 / m as f64;
        loop {
            let l = (c[seg + 1].0 - c[seg].0).hypot(c[seg + 1].1 - c[seg].1);
            if acc + l >= target || seg + 2 >= c.len() {
                let t = if l > 0.0 { (target - acc) / l } else { 0.0 };
                out.push((
                    c[seg].0 + t * (c[seg + 1].0 - c[seg].0),
                    c[seg].1 + t * (c[seg + 1].1 - c[seg].1),
                ));
                break;
            }
            acc += l;
            seg += 1;
        }
    }
    out
}

/// Revolve an (R, Z) loop (m) about the Z axis into an `n_phi`-sided torus, in mm. With `span`
/// < 2 pi the surface is open (a cut-away for pictures).
fn revolve(curve: &[(f64, f64)], n_phi: usize, span: f64, start: f64) -> IndexedMesh {
    let full = (span - 2.0 * PI).abs() < 1e-12;
    let rings = if full { n_phi } else { n_phi + 1 };
    let m = curve.len();
    let mut mesh = IndexedMesh::new();
    for k in 0..rings {
        let phi = start + span * k as f64 / n_phi as f64;
        let (s, c) = phi.sin_cos();
        for &(r, z) in curve {
            mesh.vertices
                .push(Point3::new(r * c * MM, r * s * MM, z * MM));
        }
    }
    for k in 0..n_phi {
        let k1 = if full { (k + 1) % n_phi } else { k + 1 };
        for j in 0..m {
            let j1 = (j + 1) % m;
            let (a, b, c, d) = (k * m + j, k1 * m + j, k1 * m + j1, k * m + j1);
            mesh.faces.push([a as u32, b as u32, c as u32]);
            mesh.faces.push([a as u32, c as u32, d as u32]);
        }
    }
    if full && mesh.signed_volume() < 0.0 {
        for f in &mut mesh.faces {
            f.swap(1, 2);
        }
    }
    mesh
}

/// Exact volume and area of the solid swept by revolving a polygon (Pappus), in mm^3 and mm^2.
fn pappus(curve: &[(f64, f64)]) -> (f64, f64) {
    let m = curve.len();
    let (mut moment, mut area) = (0.0, 0.0);
    for j in 0..m {
        let (a, b) = (curve[j], curve[(j + 1) % m]);
        moment += (a.0 + b.0) * (a.0 * b.1 - b.0 * a.1) / 6.0; // integral of R dA
        area += 0.5 * (a.0 + b.0) * (b.0 - a.0).hypot(b.1 - a.1); // integral of R dl
    }
    (
        2.0 * PI * moment.abs() * MM.powi(3),
        2.0 * PI * area * MM.powi(2),
    )
}

/// Mesh a solid with marching cubes and weld its coincident vertices.
fn solid_mesh(s: &Solid, cell: f64) -> IndexedMesh {
    let mut m = s.mesh(cell).geometry;
    weld_vertices(&mut m, 1e-3);
    m.faces
        .retain(|f| f[0] != f[1] && f[1] != f[2] && f[0] != f[2]);
    remove_unreferenced_vertices(&mut m);
    m
}

fn merge(parts: &[IndexedMesh]) -> IndexedMesh {
    let mut out = IndexedMesh::new();
    for p in parts {
        let base = out.vertices.len() as u32;
        out.vertices.extend_from_slice(&p.vertices);
        out.faces
            .extend(p.faces.iter().map(|f| f.map(|i| i + base)));
    }
    out
}

fn rotate_z(m: &IndexedMesh, phi: f64) -> IndexedMesh {
    let (s, c) = phi.sin_cos();
    let mut out = m.clone();
    for v in &mut out.vertices {
        *v = Point3::new(v.x * c - v.y * s, v.x * s + v.y * c, v.z);
    }
    out
}

/// A PF coil block as an axisymmetric ring solid, in mm.
fn ring(c: &Coil) -> Solid {
    let (r0, r1) = ((c.r - 0.5 * c.w) * MM, (c.r + 0.5 * c.w) * MM);
    let h = 0.5 * c.h * MM;
    let outer = Solid::cylinder(r1, h);
    let body = if r0 > 0.0 {
        outer.subtract(Solid::cylinder(r0, h + 1.0))
    } else {
        outer
    };
    body.translate(Vector3::new(0.0, 0.0, c.z * MM))
}

/// A revolve that cf-design does not ship, built on its `user_fn` escape hatch: the signed 2D
/// distance to a closed (R, Z) polygon (m), evaluated at (sqrt(x^2 + y^2), z), in mm. Exact SDF
/// wherever the profile's is (a polygon's is).
fn revolve_solid(profile: &[(f64, f64)]) -> Solid {
    let poly: Vec<(f64, f64)> = profile.iter().map(|p| (p.0 * MM, p.1 * MM)).collect();
    let rmax = poly.iter().map(|p| p.0).fold(0.0, f64::max) + 10.0;
    let (zmin, zmax) = poly
        .iter()
        .fold((f64::MAX, f64::MIN), |a, p| (a.0.min(p.1), a.1.max(p.1)));
    let bounds = cortenforge::cf_design::Aabb::new(
        Point3::new(-rmax, -rmax, zmin - 10.0),
        Point3::new(rmax, rmax, zmax + 10.0),
    );
    Solid::user_fn(
        move |p| {
            let q = (p.x.hypot(p.y), p.z);
            let n = poly.len();
            let mut d2 = f64::MAX;
            for k in 0..n {
                let (a, b) = (poly[k], poly[(k + 1) % n]);
                let (ex, ez) = (b.0 - a.0, b.1 - a.1);
                let t =
                    (((q.0 - a.0) * ex + (q.1 - a.1) * ez) / (ex * ex + ez * ez)).clamp(0.0, 1.0);
                let (dx, dz) = (q.0 - a.0 - t * ex, q.1 - a.1 - t * ez);
                d2 = d2.min(dx * dx + dz * dz);
            }
            let d = d2.sqrt();
            if inside(&poly, q) { -d } else { d }
        },
        bounds,
    )
}

/// Keep the triangles whose centroid lies outside the azimuth wedge (a0, a1): a cut-away.
fn cut(m: &IndexedMesh, a0: f64, a1: f64) -> IndexedMesh {
    let mut out = m.clone();
    out.faces.retain(|f| {
        let (mut x, mut y) = (0.0, 0.0);
        for &i in f {
            x += m.vertices[i as usize].x;
            y += m.vertices[i as usize].y;
        }
        let a = y.atan2(x);
        !(a > a0 && a < a1)
    });
    out
}

struct Report {
    rows: Vec<String>,
}

impl Report {
    fn part(&mut self, name: &str, m: &IndexedMesh, exact: Option<f64>, secs: f64) {
        let r = validate_mesh(m);
        let vol = m.signed_volume();
        let check = exact.map_or(String::from("-"), |e| format!("{:+.2e}", vol / e - 1.0));
        self.rows.push(format!(
            "  {name:<28} {:>9} faces  watertight {:<5} manifold {:<5}  volume {:>11.4e} mm^3  vs exact {check:>9}  {secs:.1} s",
            m.faces.len(),
            r.is_watertight,
            r.is_manifold,
            vol
        ));
    }
}

pub fn run() {
    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("out")
        .join("mesh");
    std::fs::create_dir_all(&out).expect("create out/mesh");
    let grid = shape::grid();
    let solver = Solver::new(&grid, true);
    let case = shape::CASES[4];
    let mut machine = shape::machine(case.1, case.2);
    let prof = shape::profile();
    let eq = solve(&grid, &solver, &machine, &prof, &Options::default(), None);
    let s = summarize(&grid, &eq, &prof);
    if let Some(i) = eq.control_current {
        machine.coils[0].current = i;
        machine.coils[1].current = i;
    }
    println!(
        "Equilibrium '{}': kappa {:.3}, delta {:+.3}, plasma volume {:.4} m^3 (node sum)",
        case.0, s.kappa, s.delta_upper, s.volume
    );
    let mut rep = Report { rows: Vec::new() };

    // Flux surfaces and the vessel wall: revolved contours.
    let n_phi = 256;
    let mut surfaces = Vec::new();
    for psin in [0.25, 0.5, 0.75, 1.0] {
        let level = eq.axis.2 - psin * (eq.axis.2 - eq.psi_b) * (1.0 - 1e-12);
        let Some(c) = contours(&grid, &eq.psi, level)
            .into_iter()
            .filter(|c| c.len() > 8 && inside(c, (eq.axis.0, eq.axis.1)))
            .min_by(|a, b| perimeter(a).total_cmp(&perimeter(b)))
        else {
            continue;
        };
        let curve = resample(&c, 256);
        let t0 = Instant::now();
        let m = revolve(&curve, n_phi, 2.0 * PI, 0.0);
        let (v, a) = pappus(&curve);
        let name = if psin < 1.0 {
            format!("flux surface psiN {psin}")
        } else {
            "last closed flux surface".into()
        };
        rep.part(&name, &m, Some(v), t0.elapsed().as_secs_f64());
        if psin == 1.0 {
            rep.rows.push(format!(
                "    area {:.4e} mm^2 vs Pappus {:.4e} ({:+.2e}); volume vs the equilibrium's node sum {:+.2e}",
                m.surface_area(),
                a,
                m.surface_area() / a - 1.0,
                m.signed_volume() / (s.volume * MM.powi(3)) - 1.0
            ));
        }
        let file = if psin < 1.0 {
            format!("surface_psin_{:.2}.stl", psin)
        } else {
            "lcfs.stl".into()
        };
        save_mesh(&m, out.join(&file)).expect("write STL");
        surfaces.push((psin, curve));
    }
    let wall: Vec<(f64, f64)> = machine.region.clone();
    let t0 = Instant::now();
    let vessel = revolve(&wall, n_phi, 2.0 * PI, 0.0);
    rep.part(
        "vessel wall (surface)",
        &vessel,
        Some(pappus(&wall).0),
        t0.elapsed().as_secs_f64(),
    );
    save_mesh(&vessel, out.join("vessel.stl")).expect("write STL");
    // The same vessel as a cf-design solid through the user_fn revolve (a finding probe).
    let t0 = Instant::now();
    let vessel_solid = solid_mesh(&revolve_solid(&wall), 10.0);
    rep.part(
        "vessel via user_fn revolve",
        &vessel_solid,
        Some(pappus(&wall).0),
        t0.elapsed().as_secs_f64(),
    );

    // PF coils and the central solenoid: cf-design rings.
    let mut pf = Vec::new();
    for c in &machine.coils {
        let t0 = Instant::now();
        let m = solid_mesh(&ring(c), 5.0);
        let (r0, r1) = ((c.r - 0.5 * c.w) * MM, (c.r + 0.5 * c.w) * MM);
        let exact = PI * (r1 * r1 - r0 * r0) * c.h * MM;
        rep.part(
            &format!("{} coil at Z {:+.2} m", c.name, c.z),
            &m,
            Some(exact),
            t0.elapsed().as_secs_f64(),
        );
        pf.push(m);
    }
    let pf_all = merge(&pf);
    save_mesh(&pf_all, out.join("pf_coils.stl")).expect("write STL");

    // TF coils: one cf-design torus in its poloidal plane, meshed once, copied round.
    let tf = TfSet::for_field(16, 1.05, 0.75, 1.0, shape::R0);
    let minor = 0.04;
    let t0 = Instant::now();
    let torus = Solid::torus(tf.radius * MM, minor * MM)
        .rotate(UnitQuaternion::from_axis_angle(
            &Vector3::x_axis(),
            0.5 * PI,
        ))
        .translate(Vector3::new(tf.rc * MM, 0.0, 0.0));
    let one = solid_mesh(&torus, 2.5);
    let exact = 2.0 * PI * PI * tf.radius * MM * (minor * MM).powi(2);
    rep.part(
        "TF coil (one of 16)",
        &one,
        Some(exact),
        t0.elapsed().as_secs_f64(),
    );
    let tf_all = merge(
        &(0..tf.n)
            .map(|k| rotate_z(&one, 2.0 * PI * k as f64 / tf.n as f64))
            .collect::<Vec<_>>(),
    );
    save_mesh(&tf_all, out.join("tf_coils.stl")).expect("write STL");

    println!(
        "Parts (mm; exact volumes: Pappus for revolved polygons, pi (r1^2 - r0^2) h for rings, 2 pi^2 R r^2 for the torus):"
    );
    for r in &rep.rows {
        println!("{r}");
    }
    println!("  STL files in {}", out.display());

    // A cut-away picture: drop the quadrant facing the camera.
    let t0 = Instant::now();
    let (a0, a1) = (-0.5 * PI, 0.0);
    let mut items = vec![
        (cut(&tf_all, a0, a1), [200u8, 120, 40]),
        (cut(&pf_all, a0, a1), [60, 110, 190]),
        (revolve(&wall, 96, 1.5 * PI, 0.0), [170, 170, 170]),
    ];
    let colors = [
        [90u8, 160, 230],
        [80, 190, 140],
        [230, 180, 60],
        [210, 60, 60],
    ];
    for (k, (_, curve)) in surfaces.iter().enumerate() {
        // Inner surfaces cut deeper so each shows.
        let span = 1.5 * PI - 0.25 * PI * (surfaces.len() - 1 - k) as f64;
        items.push((revolve(curve, 96, span, 0.0), colors[k % colors.len()]));
    }
    let scene = Scene { items };
    let path = out.join("tokamak.bmp");
    render_bmp(&scene, [-1.0, 1.0, -0.55], 1400, 1000, &path);
    println!(
        "  render: {} ({:.1} s)",
        path.display(),
        t0.elapsed().as_secs_f64()
    );
}
