#pragma once
#include "ucntrap/config.hpp"
#include "ucntrap/source/initial_condition_source.hpp"
#include <cstddef>
#include <random>

namespace ucntrap {

class RandomSource final : public InitialConditionSource {
public:
    explicit RandomSource(const SimulationConfig& config);

    bool has_next() const override;
    State next() override;

private:
    SimulationConfig config_;
    std::size_t generated_ = 0;
    mutable std::mt19937_64 rng_;
    
    double uniform(double lo, double hi);
};

} // namespace ucntrap