#pragma once
#include <random>

namespace ucntrap {
    
class RandomEngine {
public:
    explicit RandomEngine(unsigned long long seed) : rng_(seed) {}

    double uniform01() {
        return dist_(rng_);
    }

    std::mt19937_64& engine() { return rng_; }

private:
    std::mt19937_64 rng_;
    std::uniform_real_distribution<double> dist_{0.0, 1.0};
};

} // namespace ucntrap 