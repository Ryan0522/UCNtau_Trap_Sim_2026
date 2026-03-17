#pragma once
#include "ucntrap/config.hpp"
#include "ucntrap/random.hpp"

namespace ucntrap {

class Runner {
public:
    explicit Runner(SimulationConfig config);
    int run() const;

private:
    void print_config(SimulationConfig config, int size, int total_traj) const;
    SimulationConfig config_;
};

} // namespace ucntrap