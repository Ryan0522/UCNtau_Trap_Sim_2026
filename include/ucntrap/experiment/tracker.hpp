#pragma once
#include "ucntrap/state.hpp"

namespace ucntrap
{
    
struct Result {
    // time and energy
    double t_start = 0.0; // neutron entry time
    double t_final = 0.0;
    double e_start = 0.0;
    double e_final = 0.0;

    // Initial phase space coordinates
    double x_start = 0.0; double y_start = 0.0; double z_start = 0.0;
    double vx_start = 0.0; double vy_start = 0.0; double vz_start = 0.0;

    // Final phase space coordinates
    double x_final = 0.0; double y_final = 0.0; double z_final = 0.0;
    double vx_final = 0.0; double vy_final = 0.0; double vz_final = 0.0;

    // Parameter checks
    double z_off = 0.0;
    int n_hit = 0;
    int n_hit_house_low = 0;
    int n_hit_house_high = 0;
    double death_time = 0.0;

    int code = 0;
};

class Tracker {
public:
    virtual ~Tracker() = default;
    virtual Result run(const State& initial) = 0;
};

} // namespace ucntrap
