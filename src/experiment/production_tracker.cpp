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

bool ProductionTracker::check_acceptance(double x, double prev_y, double curr_y) const {
    bool in_x_range = (x > -0.3524 && x < 0.0476);

    switch (config_.dagger_mode) {
        case DaggerMode::Fast:
            return in_x_range;
        case DaggerMode::Slow:
            return in_x_range && (prev_y < 0 && curr_y > 0);
        case DaggerMode::Segmented:
            return  (x > -0.3302 && x < -0.2794) ||
                    (x > -0.2286 && x < -0.1778) ||
                    (x > -0.1270 && x < -0.0762) || 
                    (x > -0.0254 && x <  0.0254);
    }
    return false;
}

Result ProductionTracker::run(const State& startial) const {
    Result res;
    State s = startial;
    double t = 0.0;
    double energy = 0.0;

    std::vector<double> tang = {0.0, 0.0, 1.0};
    std::mt19937_64 rng(config_.seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

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
    while (t - cleaning_time < 3000.0) { // Max simulation time
        State prev_s = s;
        integrator_.step(s, t, config_.dt, field_model_);
        t += config_.dt;
        double t_exp = t - cleaning_time;

        // Check for decay
        if (t_exp > res.death_time) {
            res.code = 1; // Decay
            goto finalize;
        }

        // Check for raised cleaning
        if (std::abs(s.z - (cd::kBaseZ + config_.raised_cleaning_height)) < constants::kEpsilon) {
            res.code = -3; // Raised cleaning hit
            goto finalize;
        }

        // Check for y-plane crossing
        if ((prev_s.y < 0 && s.y > 0) || (prev_s.y > 0 && s.y < 0)) {
            // Linear inteprolation to find the exact crossing point
            double frac = std::abs(prev_s.y) / (std::abs(s.y) + std::abs(prev_s.y));
            double pred_x = prev_s.x + frac * (s.x - prev_s.x);
            double pred_z = prev_s.z + frac * (s.z - prev_s.z);
            double z_off = dagger_.get_z_offset(t_exp);

            if (dagger_.check_collision(State{pred_x, 0.0, pred_z}, t_exp)) {
                res.n_hit++;

                // 1. Acceptance Filter
                bool accepted = check_acceptance(pred_x, prev_s.y, s.y);

                if (accepted) {
                    double e_perp = (s.py * s.py) / (2.0 * constants::kMassN);
                    if (SurfaceModel::check_absorption(e_perp, config_.bthick, pred_x, pred_z, z_off)) {
                        res.code = 0; // Success: Detected
                        goto finalize;
                    } else {
                        if (dist(rng) <= config_.wall_loss_prob) {
                            res.code = -5; // Wall Loss
                            goto finalize;
                        }
                    }

                    // 4. Non-absorbing hit
                    std::vector<double> norm = {0.0, (prev_s.y > 0 ? 1.0 : -1.0), 0.0};
                    std::vector<double> tang = {0.0, 0.0, 1.0};
                    SurfaceModel::reflect(s, norm, tang);

                    // Reposition state to y = 0
                    s.x = pred_x; s.y = 0.0; s.z = pred_z;
                }
            }

            // Check for NaN
            if (std::isnan(s.x)) {
                res.code = -4;
                goto finalize;
            }
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