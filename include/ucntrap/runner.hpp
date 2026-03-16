#pragma once
#include "ucntrap/config.hpp"
#include "ucntrap/random.hpp"

namespace ucntrap {

class Runner {
public:
    explicit Runner(SimulationConfig config);
    int run() const;

private:
    SimulationConfig config_;
};

} // namespace ucntrap