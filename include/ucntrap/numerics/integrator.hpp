#pragma once
#include "ucntrap/state.hpp"
#include "ucntrap/physics/field_model.hpp"

namespace ucntrap {

class Integrator {
public:
    virtual ~Integrator() = default;
    virtual void step(State& state, double t, double dt,
                      const FieldModel& field) const = 0;

    // virtual bool step_with_defect(State& s, double t, double dt,
    //                                 const FieldModel& field, double defect_prob,
    //                                 RandomEngine& rng, double& total_energy) const = 0;
};

const Integrator& default_integrator();

} // namespace ucntrap