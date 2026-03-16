#include "ucntrap/experiment/production_tracker.hpp"
#include "ucntrap/physics/surface_model.hpp"
#include "ucntrap/constants.hpp"

#include <cmath>
#include <limits>
#include <random>

namespace ucntrap {

namespace cd = constants::dagger;

ProductionTracker::ProductionTracker(const SimulationConfig& config,
                                     const FieldModel& field_model,
                                     const Integrator& integrator,
                                     const Dagger& dagger)
    : config_(config),
      field_model_(field_model),
      integrator_(integrator),
      dagger_(dagger) {}

Result ProductionTracker::run(const State& startial) const {
    Result res;
    State s = startial;
    double t = 0.0;
    double energy = 0.0;

    // 1. Record startial conditions
    res.t_start = 0.0;
    res.x_start = s.x; res.y_start = s.y; res.z_start = s.z;
    res.vx_start = s.px / constants::kMassN; 
    res.vy_start = s.py / constants::kMassN; 
    res.vz_start = s.pz / constants::kMassN;

    double u_start = field_model_.potential(s, t);
    res.e_start = u_start - constants::kMinU + (s.px*s.px + s.py*s.py + s.pz*s.pz) / (2.0 * constants::kMassN);

    // Set death time and clean time
    res.death_time = 9999;
    double cleaning_time = config_.cleaning_time;

    // 2. Cleaning Phase
    int cleaning_steps = static_cast<int>(cleaning_time / config_.dt);
    for (int i = 0; i < cleaning_steps; ++i) {
        State prev_s = s;
        integrator_.step(s, t, config_.dt, field_model_);
        t += config_.dt;

        // Check if cleaned
        if ((prev_s.z < (cd::kBaseZ + config_.cleaning_height) && s.z > (cd::kBaseZ + config_.cleaning_height)) ||
            (prev_s.z > (cd::kBaseZ + config_.cleaning_height) && s.z < (cd::kBaseZ + config_.cleaning_height))) {
            if (s.py > 0) {
                res.code = -2; // Cleaned
                goto finalize;
            }
        }
    }
    
    // 3. Detection Loop
    while (true) {
        State prev_s = s;
        integrator_.step(s, t, config_.dt, field_model_);
        t += config_.dt;

        double current_t_exp = t - cleaning_time;

        // Check for decay
        if (current_t_exp > res.death_time) {
            res.code = 1; // Decay
            goto finalize;
        }

        // Check for raised cleaning
        if (std::abs(s.z - (cd::kBaseZ + config_.raised_cleaning_height)) < constants::kEpsilon) {
            res.code = -3; // Raised cleaning hit
            goto finalize;
        }

        // Check for dagger hit
        if (dagger_.check_collision(s, current_t_exp)) {
            res.n_hit++;
            res.code = 1; // Hit
            goto finalize;
        }

        // Check if Nan
        if (std::isnan(s.x)) {
            res.code = -4;
            goto finalize;
        }
    }

finalize:
    // 4. Record final conditions
    res.t_final = t - cleaning_time;
    res.x_final = s.x; res.y_final = s.y; res.z_final = s.z;
    res.vx_final = s.px / constants::kMassN;
    res.vy_final = s.py / constants::kMassN;
    res.vz_final = s.pz / constants::kMassN;
    res.e_final = field_model_.potential(s, t) - constants::kMinU + (s.px*s.px + s.py*s.py + s.pz*s.pz) / (2.0 * constants::kMassN);
    res.z_off = dagger_.get_z_offset(res.t_final);

    return res;
}   

} // namespace ucntrap