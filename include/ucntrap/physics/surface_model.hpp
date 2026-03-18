#pragma once
#include "ucntrap/state.hpp"
#include "ucntrap/random.hpp"
#include <array>

namespace ucntrap {

class SurfaceModel {
public:
    using Vec3 = std::array<double, 3>;

    static bool check_absorption(double e_perp, double b_thick, double x, double z, double z_off, double zeta_cut, RandomEngine& rng);

    static void reflect(State& s, const Vec3& norm, const Vec3& tang, RandomEngine& rng);

private:
    static double calculate_absorption_prob(double e_perp, double b_thick);
};

} // namespace ucntrap