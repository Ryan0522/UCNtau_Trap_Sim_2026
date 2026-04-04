#include "ucntrap/experiment/energy_tracker.hpp"
#include "ucntrap/physics/surface_model.hpp"
#include "ucntrap/constants.hpp"

#include <cmath>
#include <limits>

namespace ucntrap {

Result EnergyTracker::run(const State& initial) {
    Result res;
    State s = initial;
    double t = 0.0;

    // 1. Record initial conditions
    res.t_start = 0.0;
    res.x_start = s.x; res.y_start = s.y; res.z_start = s.z;
    res.vx_start = s.px / constants::kMassN; 
    res.vy_start = s.py / constants::kMassN; 
    res.vz_start = s.pz / constants::kMassN;

    double u_start = field_model_.potential(s, t);
    res.e_start = u_start - constants::kMinU + (s.px*s.px + s.py*s.py + s.pz*s.pz) / (2.0 * constants::kMassN);
    
    // 2. Simulate until death or timeout
    while (t < config_.hold_time) {
        integrator_.step(s, t, config_.dt, field_model_);
        t += config_.dt;
    }

    res.t_final = t;
    res.x_final = s.x; res.y_final = s.y; res.z_final = s.z;
    res.vx_final = s.px / constants::kMassN;
    res.vy_final = s.py / constants::kMassN;
    res.vz_final = s.pz / constants::kMassN;
    res.e_final = field_model_.potential(s, t) - constants::kMinU + (s.px*s.px + s.py*s.py + s.pz*s.pz) / (2.0 * constants::kMassN);
    res.code = 0;

    return res;
}

} // namespace ucntrap