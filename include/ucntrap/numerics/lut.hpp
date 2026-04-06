#pragma once
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

    // for sin/cos periodic table assume n = 65536 and input is in [0, 2pi)
    inline double eval_periodic(double x) const {
        double pos = (x - min_) * inv_step_;
        
        int64_t i0 = static_cast<int64_t>(std::floor(pos));
        double frac = pos - static_cast<double>(i0);

        uint64_t mask = n_ - 1;
        uint64_t idx0 = static_cast<uint64_t>(i0) & mask;
        uint64_t idx1 = (idx0 + 1) & mask;

        // 4. Linear interpolation: y = y0 + frac * (y1 - y0)
        return table_[idx0] + frac * (table_[idx1] - table_[idx0]);
    }

private:
    std::vector<double> table_;
    double min_ = 0, max_ = 0, step_ = 0, inv_step_ = 0;
    std::size_t n_ = 0;
};

} // namespace ucntrap