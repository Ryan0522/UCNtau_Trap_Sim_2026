#pragma once
#include "ucntrap/experiment/tracker.hpp"
#include "ucntrap/physics/field_model.hpp"
#include "ucntrap/numerics/integrator.hpp"
#include "ucntrap/config.hpp"

namespace ucntrap {

class EnergyTracker : public Tracker {
public:
    EnergyTracker(const SimulationConfig& config,
                      const FieldModel& field_model,
                      const Integrator& integrator)
        : config_(config), field_model_(field_model), integrator_(integrator) {}

    Result run(const State& initial) override;

private:
    SimulationConfig config_;
    const FieldModel& field_model_;
    const Integrator& integrator_;
};

} // namespace ucntrap