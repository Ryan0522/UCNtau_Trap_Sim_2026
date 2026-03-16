#include "ucntrap/physics/dagger.hpp"
#include "ucntrap/constants.hpp"

#include <cmath>
#include <algorithm>

namespace ucntrap {

namespace cd = constants::dagger;

Dagger::Dagger(const std::vector<double>& dip_heights, const std::vector<double>& dip_ends)
    : dip_heights_(dip_heights), dip_ends_(dip_ends) {}

double Dagger::get_z_offset(double t) const {
    if (dip_ends_.empty() || dip_heights_.empty()) {
        return 0.0; // No dips defined
    }

    if (t > dip_ends_.back()) {
        return 0.01;
    }

    std::size_t i = 0;
    for (; i < dip_ends_.size(); ++i) {
        if (dip_ends_[i] > t) break;
    }

    if (i == 0) return dip_heights_[0];

    double target = dip_heights_[i];
    double start = dip_heights_[i - 1];
    double distance_um = std::abs(target - start) * 1e6;

    double move_time = t - dip_ends_[i - 1];
    double sign = (target - start > 0) ? 1.0 : (target - start < 0 ? -1.0 : 0.0);

    double plateau_time = (distance_um - cd::kVel * (cd::kVel / cd::kAcc) * cd::kSpr) / (cd::kVel * cd::kSpr);

    if (plateau_time > 0) {
       if (move_time < (cd::kVel / cd::kAcc)) {
            return start + sign * (0.5 * cd::kAcc * cd::kSpr * move_time * move_time) / 1000000.0;
        } 
        else if (move_time < (cd::kVel / cd::kAcc + plateau_time)) {
            return start + sign * (0.5 * (cd::kVel / cd::kAcc) * cd::kVel * cd::kSpr + 
                                  (move_time - cd::kVel / cd::kAcc) * cd::kVel * cd::kSpr) / 1000000.0;
        } 
        else if (move_time < (2.0 * (cd::kVel / cd::kAcc) + plateau_time)) {
            double dt_dec = move_time - plateau_time - (cd::kVel / cd::kAcc);
            return start + sign * (0.5 * (cd::kVel / cd::kAcc) * cd::kVel * cd::kSpr + plateau_time * cd::kVel * cd::kSpr +
                                  cd::kVel * cd::kSpr * dt_dec - 0.5 * cd::kAcc * cd::kSpr * dt_dec * dt_dec) / 1000000.0;
        }
    } else {
        double half_time = std::sqrt(distance_um / (cd::kAcc * cd::kSpr));
        if (move_time < half_time) {
            return start + sign * (0.5 * cd::kAcc * cd::kSpr * move_time * move_time) / 1000000.0;
        } else if (move_time < 2.0 * half_time) {
            double dt_dec = move_time - half_time;
            return start + (target - start) / 2.0 + 
                   sign * (cd::kAcc * cd::kSpr * half_time * dt_dec - 0.5 * cd::kAcc * cd::kSpr * dt_dec * dt_dec) / 1000000.0;
        }
    }
    
    return target;
}

bool Dagger::check_collision(const State& s, double t) const {
    const double z_off = get_z_offset(t);
    const double x = s.x;
    const double z = s.z;

    // 1. Calculate zeta (distance from neutron to dagger surface)
    double zeta = 0.0;
    const double z_rel = std::abs(z - z_off);
    if (x > 0.0) {
        zeta = 0.5 - std::sqrt(x * x + std::pow(z_rel - 1.0, 2));
    } else {
        zeta = 1.0 - std::sqrt(x * x + std::pow(z_rel - 0.5, 2));
    }

    // 2. Check if hit is on Dagger
    if (x > cd::kXMin && x < cd::kXMax && zeta > 0.0 && z < (cd::kBaseZ + z_off + cd::kDaggerZShift)) {
        return true;
    }

    // 3. Check if hit is on House
    const double z_house_base = cd::kBaseZ + z_off + cd::kDaggerZShift;
    const double dx = std::abs(x - cd::kHouseXCenter);

    // House Hit Low
    if (z >= z_house_base && z < (z_house_base + cd::kHouseLowH)) {
        if (dx < (0.40 + 2.0179 * (z - z_house_base)) / 2.0) return true;
    }
    // House Hit High
    if (z >= (z_house_base + cd::kHouseLowH) && z < (z_house_base + cd::kHouseHighH)) {
        if (dx < 0.69215 / 2.0) return true;
    }

    return false;
}

HitInfo Dagger::classify_crossing(double x, double z, double t_exp) const {
    double z_off = get_z_offset(t_exp);

    if (x > cd::kXMin && x < cd::kXMax) {
        double z_rel = std::abs(z - z_off);
        double zeta = (x > 0) ? (0.5 - std::sqrt(x*x + std::pow(z_rel - 1.0, 2))) 
                              : (1.0 - std::sqrt(x*x + std::pow(z_rel - 0.5, 2)));
        
        if (zeta > 0.0 && z < (cd::kBaseZ + z_off + cd::kDaggerZShift)) {
            return {HitType::Dagger, x, z, z_off};
        }
    }

    if (check_house_hit_low(x, z, z_off)) return {HitType::HouseLow, x, z, z_off};
    if (check_house_hit_high(x, z, z_off)) return {HitType::HouseHigh, x, z, z_off};

    return {HitType::None};
}

bool Dagger::check_house_hit_low(double x, double z, double z_off) const {
    const double z_house_base = cd::kBaseZ + z_off + cd::kDaggerZShift;
    const double dx = std::abs(x - cd::kHouseXCenter);

    if (z >= z_house_base && z < (z_house_base + cd::kHouseLowH)) {
        return dx < (0.40 + 2.0179 * (z - z_house_base)) / 2.0;
    }
    return false;
}

bool Dagger::check_house_hit_high(double x, double z, double z_off) const {
    const double z_house_base = cd::kBaseZ + z_off + cd::kDaggerZShift;
    const double dx = std::abs(x - cd::kHouseXCenter);

    if (z >= (z_house_base + cd::kHouseLowH) &&
        z <  (z_house_base + cd::kHouseHighH)) {
        return dx < 0.69215 / 2.0;
    }
    return false;
}

} // namespace ucntrap