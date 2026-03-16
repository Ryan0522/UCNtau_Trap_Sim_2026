#pragma once

#include "ucntrap/state.hpp"

namespace ucntrap
{
    
struct Result {
    double death_time = 0.0;
    int code = 0;
};

class Tracker {
public:
    virtual ~Tracker() = default;
    virtual Result run(const State& initial) const = 0;
};

} // namespace ucntrap
