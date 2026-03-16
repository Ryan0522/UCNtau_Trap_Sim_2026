#pragma once

#include "ucntrap/state.hpp"
#include "ucntrap/physics/field_model.hpp"

namespace ucntrap {

class Integrator {
public:
    virtual ~Integrator() = default;
    virtual void step(State& state, double t, double dt,
                        const FieldModel& field) const = 0;
};

} // namespace ucntrap