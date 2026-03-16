#pragma once
#include "ucntrap/experiment/tracker.hpp"
#include "ucntrap/physics/dagger.hpp"
#include "ucntrap/physics/field_model.hpp"
#include "ucntrap/numerics/integrator.hpp"
#include "ucntrap/config.hpp"

namespace ucntrap {

class ProductionTracker : public Tracker {
public:
    ProductionTracker(const SimulationConfig& config,
                      const FieldModel& field_model,
                      const Integrator& integrator,
                      const Dagger& dagger);

    Result run(const State& initial) const override;

private:
    SimulationConfig config_;
    const FieldModel& field_model_;
    const Integrator& integrator_;
    const Dagger& dagger_;
};

} // namespace ucntrap