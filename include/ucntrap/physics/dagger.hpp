#pragma once
#include "ucntrap/state.hpp"
#include <vector>

namespace ucntrap {

class Dagger {
public:
    Dagger(const std::vector<double>& dip_heights, const std::vector<double>& dip_ends);

    // Get the z-offset for a given time t based on the defined dips
    double get_z_offset(double t) const;

    // Check if the dagger is hit by neutron (both Dagger and House)
    bool check_collision(const State& s, double t) const;

private:
    std::vector<double> dip_heights_;
    std::vector<double> dip_ends_;
};

}