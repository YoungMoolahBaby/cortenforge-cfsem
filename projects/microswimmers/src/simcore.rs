//! Probes of CortenForge's own fluid model (sim-core, MuJoCo's two models): a free sphere pushed
//! and twisted, and a free helix of capsules twisted about its axis.

use nalgebra::{Vector3, Vector6};

use cortenforge::sim::core::Model;
use cortenforge::sim::mjcf::load_model;

use crate::helix::{Helix, V3};

/// Steady state of a free body under a constant applied torque and force (world frame): its
/// linear and angular velocity (world frame) and the steps taken.
pub struct Steady {
    pub v: V3,
    pub w: V3,
    pub steps: usize,
}

/// Step until the motion stops changing (relative 1e-10 over 1000 steps), the spin scaled by the
/// body size `ell` so a body that only turns settles too; at most 50,000 steps (a free helix
/// wobbles as it turns, so its spin never settles to 1e-10: it is read after the cap).
fn settle(model: &Model, torque: V3, force: V3, ell: f64) -> Steady {
    let mut data = model.make_data();
    // sim-core's layout: [torque, force] (forward/acceleration.rs:131-132).
    data.xfrc_applied[1] = Vector6::new(
        torque.x, torque.y, torque.z, force.x, force.y, force.z,
    );
    let read = |d: &cortenforge::sim::core::Data| {
        let v = Vector3::new(d.qvel[0], d.qvel[1], d.qvel[2]);
        // A free joint's angular velocity is in the body frame.
        let w_local = Vector3::new(d.qvel[3], d.qvel[4], d.qvel[5]);
        (v, d.xquat[1] * w_local)
    };
    let mut last = (V3::zeros(), V3::zeros());
    let mut steps = 0;
    loop {
        for _ in 0..1000 {
            data.step(model).expect("step");
        }
        steps += 1000;
        let now = read(&data);
        let change = ((now.0 - last.0).norm() + ell * (now.1 - last.1).norm())
            / (now.0.norm() + ell * now.1.norm()).max(1e-300);
        last = now;
        if (change < 1e-10 && steps >= 5000) || steps >= 50_000 {
            break;
        }
    }
    Steady {
        v: last.0,
        w: last.1,
        steps,
    }
}

fn option(dt: f64, mu: f64, rho_fluid: f64) -> String {
    format!(
        r#"<option timestep="{dt}" gravity="0 0 0" viscosity="{mu}" density="{rho_fluid}"/>"#
    )
}

/// A free sphere of radius `a` in fluid of viscosity `mu`; `shape` "ellipsoid" (per-geom model)
/// or "none" (the body-level inertia-box model). Returns the steady states under a force along x
/// and a torque about z.
pub fn sphere(a: f64, mu: f64, shape: &str, dt: f64, force: f64, torque: f64) -> (Steady, Steady) {
    let xml = format!(
        r#"<mujoco model="sphere">
  {}
  <worldbody>
    <body name="ball" pos="0 0 0">
      <joint type="free"/>
      <geom type="sphere" size="{a}" density="1000" fluidshape="{shape}" contype="0" conaffinity="0"/>
    </body>
  </worldbody>
</mujoco>"#,
        option(dt, mu, 0.0)
    );
    let model = load_model(&xml).expect("sphere model loads");
    (
        settle(&model, V3::zeros(), V3::new(force, 0.0, 0.0), a),
        settle(&model, V3::new(0.0, 0.0, torque), V3::zeros(), a),
    )
}

/// A free helix of `per_turn` capsules a turn, twisted about z by `torque`.
pub fn helix(
    h: &Helix,
    per_turn: usize,
    mu: f64,
    rho_fluid: f64,
    shape: &str,
    dt: f64,
    torque: f64,
) -> Steady {
    let n = (per_turn as f64 * h.turns).round() as usize;
    let c = h.centre();
    let mut geoms = String::new();
    for i in 0..n {
        let p = h.point(h.t_max() * i as f64 / n as f64) - c;
        let q = h.point(h.t_max() * (i + 1) as f64 / n as f64) - c;
        geoms += &format!(
            "      <geom type=\"capsule\" fromto=\"{} {} {} {} {} {}\" size=\"{}\" density=\"1250\" fluidshape=\"{shape}\" contype=\"0\" conaffinity=\"0\"/>\n",
            p.x, p.y, p.z, q.x, q.y, q.z, h.wire
        );
    }
    let xml = format!(
        r#"<mujoco model="helix">
  {}
  <worldbody>
    <body name="helix" pos="0 0 0">
      <joint type="free"/>
{geoms}    </body>
  </worldbody>
</mujoco>"#,
        option(dt, mu, rho_fluid)
    );
    let model = load_model(&xml).expect("helix model loads");
    settle(&model, V3::new(0.0, 0.0, torque), V3::zeros(), h.radius)
}
