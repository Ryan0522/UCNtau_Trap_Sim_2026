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

} // namespace ucntrap
