#pragma once
#include "ucntrap/state.hpp"
#include <vector>

namespace ucntrap {

class SurfaceModel {
public:
    static bool check_absorption(double e_perp, double b_thick, double x, double z, double z_off);

    static void reflect(State& s, const std::vector<double>& norm, const std::vector<double>& tang);

private:
    static double calculate_absorption_prob(double e_perp, double b_thick);
};

} // namespace ucntrap