#pragma once
#include "ucntrap/physics/field_model.hpp"
#include <vector>

namespace ucntrap {

class TrapHalbachField final : public FieldModel {
public:
    TrapHalbachField(double heat_mult,
                     const std::vector<double>& tx,
                     const std::vector<double>& ty,
                     const std::vector<double>& tz);
    
    Force force(const State& s, double t) const override;
    double potential(const State& s, double t) const override;

private:
    double heat_mult_;
    
    const std::vector<double>& tx_, ty_, tz_;

    void get_shifted_coords(const State& s, double t, double& x, double& y, double& z) const;
};

} // namespace ucntrap