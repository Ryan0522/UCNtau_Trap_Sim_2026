#pragma once

#include "ucntrap/state.hpp"

namespace ucntrap
{

struct Force
{
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
};

class FieldModel {
public:
    virtual ~FieldModel() = default;
    virtual Force force(const State& s, double t) const = 0;
    virtual double potential(const State& s, double t) const = 0;
};

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
