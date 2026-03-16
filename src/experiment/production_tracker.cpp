#include "ucntrap/experiment/production_tracker.hpp"
#include "ucntrap/physics/surface_model.hpp"
#include "ucntrap/constants.hpp"

#include <cmath>
#include <limits>

namespace ucntrap {

namespace cd = constants::dagger;

ProductionTracker::ProductionTracker(const SimulationConfig& config,
                                     const FieldModel& field_model,
                                     const Integrator& integrator,
                                     const Dagger& dagger,
                                     RandomEngine& rng)
    : config_(config),
      field_model_(field_model),
      integrator_(integrator),
      dagger_(dagger),
      rng_(rng) {}

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

Result ProductionTracker::run(const State& initial) const {
    Result res;
    State s = initial;
    double t = 0.0;
    double energy = 0.0;

    // 1. Record initial conditions
    res.t_start = 0.0;
    res.x_start = s.x; res.y_start = s.y; res.z_start = s.z;
    res.vx_start = s.px / constants::kMassN; 
    res.vy_start = s.py / constants::kMassN; 
    res.vz_start = s.pz / constants::kMassN;

    double u_start = field_model_.potential(s, t);
    res.e_start = u_start - constants::kMinU + (s.px*s.px + s.py*s.py + s.pz*s.pz) / (2.0 * constants::kMassN);

    // Set death time and clean time
    res.death_time = 9999;
    const double cleaning_time = config_.cleaning_time;

    State prev_s = s;

    // 2. Cleaning Phase
    const int cleaning_steps = static_cast<int>(cleaning_time / config_.dt);
    for (int i = 0; i < cleaning_steps; ++i) {
        prev_s = s;
        integrator_.step(s, t, config_.dt, field_model_);
        t += config_.dt;

        // Check if cleaned
        const double clean_z = cd::kBaseZ + config_.cleaning_height;
        const bool crossed_clean =
            (prev_s.z < clean_z && s.z > clean_z) ||
            (prev_s.z > clean_z && s.z < clean_z);

        if (crossed_clean && s.py > 0.0) {
            res.code = -2;
            goto finalize;
        }

        if (std::isnan(s.x) || std::isnan(s.y) || std::isnan(s.z)) {
            res.code = -4;
            goto finalize;
        }
    }
    
    // 3. Detection Loop
    while (t - cleaning_time < 3000.0) { // Max simulation time
        prev_s = s;
        integrator_.step(s, t, config_.dt, field_model_);
        t += config_.dt;
        
        const double t_exp = t - cleaning_time;

        // Check for decay
        if (t_exp > res.death_time) {
            res.code = 1; // Decay
            goto finalize;
        }

        // Check for raised cleaning
        const double raised_clean_z = cd::kBaseZ + config_.raised_cleaning_height;
        const bool crossed_raised_clean =
            (prev_s.z < raised_clean_z && s.z > raised_clean_z) ||
            (prev_s.z > raised_clean_z && s.z < raised_clean_z);

        if (crossed_raised_clean) {
            res.code = -3; // raised clean
            goto finalize;
        }

        // Check for xz-plane crossing
        const bool crossed_y0 =
            (prev_s.y < 0.0 && s.y > 0.0) ||
            (prev_s.y > 0.0 && s.y < 0.0);
        
        if (crossed_y0) {
            // Linear inteprolation to find the exact crossing point
            const double frac = std::abs(prev_s.y) / (std::abs(s.y) + std::abs(prev_s.y));
            const double pred_x = prev_s.x + frac * (s.x - prev_s.x);
            const double pred_z = prev_s.z + frac * (s.z - prev_s.z);
            
            const HitInfo hit = dagger_.classify_crossing(pred_x, pred_z, t_exp);

            if (hit.type != HitType::None) {
                bool terminate = false;

                if (hit.type == HitType::Dagger) {
                    res.n_hit++;

                    const bool accepted = check_acceptance(hit.x, prev_s.y, s.y);

                    if (accepted) {
                        const double e_perp =
                            (prev_s.py * prev_s.py) / (2.0 * constants::kMassN);

                        if (SurfaceModel::check_absorption(
                                e_perp, config_.bthick, hit.x, hit.z, hit.z_off, rng_)) {
                            s = prev_s;
                            s.x = hit.x;
                            s.y = 0.0;
                            s.z = hit.z;
                            res.code = 0; // detected
                            terminate = true;
                        } else if (rng_.uniform01() <= config_.wall_loss_prob) {
                            s = prev_s;
                            s.x = hit.x;
                            s.y = 0.0;
                            s.z = hit.z;
                            res.code = -5; // wall loss after accepted dagger hit
                            terminate = true;
                        }
                    } else {
                        if (rng_.uniform01() <= config_.wall_loss_prob) {
                            s = prev_s;
                            s.x = hit.x;
                            s.y = 0.0;
                            s.z = hit.z;
                            res.code = -5; // rejected geometry hit lost on wall
                            terminate = true;
                        }
                    }
                } else if (hit.type == HitType::HouseLow) {
                    res.n_hit_house_low++;

                    if (rng_.uniform01() <= config_.wall_loss_prob) {
                        s = prev_s;
                        s.x = hit.x;
                        s.y = 0.0;
                        s.z = hit.z;
                        res.code = -6;
                        terminate = true;
                    }
                } else if (hit.type == HitType::HouseHigh) {
                    res.n_hit_house_high++;

                    if (rng_.uniform01() <= config_.wall_loss_prob) {
                        s = prev_s;
                        s.x = hit.x;
                        s.y = 0.0;
                        s.z = hit.z;
                        res.code = -7;
                        terminate = true;
                    }
                }

                if (terminate) {
                    goto finalize;
                }

                // Surviving physical hit: reflect from the pre-hit state.
                s = prev_s;

                const std::vector<double> norm =
                    (prev_s.y > 0.0 && prev_s.py < 0.0)
                        ? std::vector<double>{0.0, 1.0, 0.0}
                        : std::vector<double>{0.0, -1.0, 0.0};

                const std::vector<double> tang = {0.0, 0.0, 1.0};

                SurfaceModel::reflect(s, norm, tang, rng_);

                // Put state back onto the crossing plane.
                s.x = hit.x;
                s.y = 0.0;
                s.z = hit.z;
            }
        }

        if (std::isnan(s.x) || std::isnan(s.y) || std::isnan(s.z) ||
            std::isnan(s.px) || std::isnan(s.py) || std::isnan(s.pz)) {
            res.code = -4;
            goto finalize;
        }
    }

    // timeout / max sim time
    res.code = 2;

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