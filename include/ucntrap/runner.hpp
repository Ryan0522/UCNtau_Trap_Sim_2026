#pragma once
#include "ucntrap/config.hpp"

namespace ucntrap {

class Runner {
public:
    explicit Runner(SimulationConfig config);
    int run() const;

private:
    SimulationConfig config_;
};

} // namespace ucntrap