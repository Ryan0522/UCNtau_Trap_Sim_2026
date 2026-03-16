#include "ucntrap/source/random_source.hpp"
#include "ucntrap/constants.hpp"

#include <cmath>
#include <stdexcept>

namespace ucntrap {

RandomSource::RandomSource(const SimulationConfig& config)
    : config_(config), rng_(config.seed) {}

bool RandomSource::has_next() const {
    return generated_ < config_.ntraj;
}

double RandomSource::uniform(double lo, double hi) {
    std::uniform_real_distribution<double> dist(lo, hi);
    return dist(rng_);
}

State RandomSource::next() {
    if (!has_next()) {
        throw std::out_of_range("RandomSource::next() called after exhaustion");
    }

    ++generated_;

    State s{};

    s.x = uniform(-0.02, 0.02);
    s.y = uniform(-0.02, 0.02);
    s.z = uniform(0.00, 0.10);

    s.px = uniform(-1.0, 1.0);
    s.py = uniform(-1.0, 1.0);
    s.pz = uniform(-1.0, 1.0);

    return s;
}

} // namespace ucntrap