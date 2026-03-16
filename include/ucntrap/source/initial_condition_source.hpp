#pragma once

#include "ucntrap/state.hpp"

namespace ucntrap {

class InitialConditionSource {
public:
    virtual ~InitialConditionSource() = default;
    virtual bool has_next() const = 0;
    virtual State next() = 0;
};

} // namespace ucntrap