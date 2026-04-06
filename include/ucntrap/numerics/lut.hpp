#pragma once
#include "ucntrap/constants.hpp"
#include <vector>
#include <cmath>
#include <functional>
#include <cstdint>

namespace ucntrap {

class LinearLUT {
public:
    LinearLUT() = default;

    void initialize(double min, double max, std::size_t n, const std::function<double(double)>& func) {
        min_ = min;
        max_ = max;
        n_ = n;
        step_ = (max - min) / static_cast<double>(n - 1);
        inv_step_ = 1.0 / step_;

        table_.clear();
        table_.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            table_.push_back(func(min + i * step_));
        }
    }

    inline double eval(double x) const {
        if (x <= min_) return table_.front();
        if (x >= max_) return table_.back();
        double pos = (x - min_) * inv_step_;
        std::size_t i = static_cast<std::size_t>(pos);
        double frac = pos - static_cast<double>(i);

        // Linear interpolation: y = y0 + frac * (y1 - y0)
        return table_[i] + frac * (table_[i + 1] - table_[i]);
    }

    // Assuming n = 65536, min = 0, max = 2*PI
    double scale = 65536.0 / (2.0 * constants::kPi);
    uint64_t mask = 65535; // n - 1

    // for sin/cos periodic table
    inline double eval_periodic(double x) const {
        uint64_t idx = static_cast<uint64_t>(x * scale) & mask;
        return table_[idx];
    }

private:
    std::vector<double> table_;
    double min_ = 0, max_ = 0, step_ = 0, inv_step_ = 0;
    std::size_t n_ = 0;
};

} // namespace ucntrap