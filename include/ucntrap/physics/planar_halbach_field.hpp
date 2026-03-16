#pragma once

#include "ucntrap/physics/field_model.hpp"

namespace ucntrap {

class PlanarHalbachField final : public FieldModel {
public:
    PlanarHalbachField(double b_rem, double mag_space, double mag_thick, int n_terms);
    
    Force force(const State& s, double t) const override;
    double potential(const State& s, double t) const override;

private:
    double b_rem_;
    double mag_space_;
    double mag_thick_;
    int n_terms_;
};

} // namespace ucntrap