#pragma once
#include "ucntrap/state.hpp"
#include <vector>

namespace ucntrap {

enum class HitType { None, Dagger, HouseLow, HouseHigh };

struct HitInfo {
    HitType type = HitType::None;
    double x = 0.0;
    double z = 0.0;
    double z_off = 0.0;
};

class Dagger {
public:
    Dagger(const std::vector<double>& dip_heights, const std::vector<double>& dip_ends);

    // Get the z-offset for a given time t based on the defined dips
    double get_z_offset(double t) const;

    // Check if the dagger is hit by neutron (both Dagger and House)
    bool check_collision(const State& s, double t) const;

    HitInfo classify_crossing(double x, double z, double t_exp) const;

private:
    bool check_house_hit_low(double x, double z, double z_off) const;
    bool check_house_hit_high(double x, double z, double z_off) const;
    
    std::vector<double> dip_heights_;
    std::vector<double> dip_ends_;
};

}