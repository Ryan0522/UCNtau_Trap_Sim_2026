#pragma once
#include <cstddef>
#include <vector>

namespace ucntrap {

struct Trace {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    [[nodiscard]] bool empty() const noexcept {
        return x.empty() || y.empty() || z.empty();
    }

    [[nodiscard]] std::size_t size() const noexcept {
        return x.size();
    }
};

} // namespace ucntrap