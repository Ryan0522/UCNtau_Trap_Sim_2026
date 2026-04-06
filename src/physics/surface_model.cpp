#include "ucntrap/physics/surface_model.hpp"
#include "ucntrap/constants.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <vector>

namespace ucntrap {

namespace cs = constants::surface;
using Vec3 = std::array<double, 3>;

namespace {
using ComplexMat = std::array<std::complex<double>, 4>;

ComplexMat matmul(const ComplexMat& a, const ComplexMat& b) {
    return {
        a[0] * b[0] + a[1] * b[2],
        a[0] * b[1] + a[1] * b[3],
        a[2] * b[0] + a[3] * b[2],
        a[2] * b[1] + a[3] * b[3]
    };
}

std::complex<double> k_wave(double e_perp, std::complex<double> u) {
    return std::sqrt(
        (2.0 * constants::kMassN / (constants::kHbar * constants::kHbar)) *
        (e_perp - u)
    );
}

ComplexMat m_matrix(std::complex<double> kn,
                    std::complex<double> kn_minus_1,
                    double z) {
    const std::complex<double> gamma = kn_minus_1 / kn;
    const std::complex<double> i(0.0, 1.0);

    ComplexMat res;
    res[0] = 0.5 * (1.0 + gamma) * std::exp(i * (kn_minus_1 - kn) * z);
    res[1] = 0.5 * (1.0 - gamma) * std::exp(-i * (kn_minus_1 + kn) * z);
    res[2] = 0.5 * (1.0 - gamma) * std::exp(i * (kn_minus_1 + kn) * z);
    res[3] = 0.5 * (1.0 + gamma) * std::exp(-i * (kn_minus_1 - kn) * z);
    return res;
}

std::vector<double> cross3(const Vec3& a, const Vec3& b) {
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

double clamp01(double x) {
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}
} // namespace

double SurfaceModel::calculate_absorption_prob(double e_perp, double b_thick) {
    static const double v_boron =
        (2.0 * constants::kPi * (constants::kHbar * constants::kHbar) /
         constants::kMassN) *
        cs::kABoron * cs::kNBoron;

    static const double w_boron =
        (constants::kHbar / 2.0) * cs::kNBoron * cs::kSigmaBoron;

    static const double v_zns =
        (2.0 * constants::kPi * (constants::kHbar * constants::kHbar) /
         constants::kMassN) *
            cs::kAZinc * cs::kNZinc +
        (2.0 * constants::kPi * (constants::kHbar * constants::kHbar) /
         constants::kMassN) *
            cs::kASulfur * cs::kNSulfur;

    static const double w_zns =
        (constants::kHbar / 2.0) *
        (cs::kNZinc * cs::kSigmaZinc + cs::kNSulfur * cs::kSigmaSulfur);

    static const std::vector<std::complex<double>> potentials = {
        {0.0, 0.0},
        {v_boron, -w_boron},
        {v_zns, -w_zns}
    };

    // b_thick is in nm
    const std::vector<double> z_boundaries = {
        0.0,
        b_thick * 1.0e-9,
        10000.0e-9
    };

    ComplexMat m_bar = {1.0, 0.0, 0.0, 1.0};

    for (int i = static_cast<int>(potentials.size()) - 1; i > 0; --i) {
        const auto kn = k_wave(e_perp, potentials[i]);
        const auto kn_m1 = k_wave(e_perp, potentials[i - 1]);
        m_bar = matmul(m_bar, m_matrix(kn, kn_m1, z_boundaries[i - 1]));
    }

    const std::complex<double> r = -m_bar[2] / m_bar[3];
    return clamp01(1.0 - std::norm(r));
}

bool SurfaceModel::check_absorption(double e_perp,
                                    double b_thick,
                                    double x,
                                    double z,
                                    double z_off,
                                    double zeta_cut,
                                    RandomEngine& rng) {
    const double z_rel = std::abs(z - z_off);

    const double dz = (x > 0.0) ? (z_rel - 1.0) : (z_rel - 0.5);
    const double r = (x > 0.0) ? 0.5 : 1.0;
    double zeta = r - std::sqrt(x * x + dz * dz);

    double coverage = 0.0;
    if (zeta >= zeta_cut) {
        coverage = 1.0;
    } else if (zeta > 0.0) {
        coverage = zeta / zeta_cut;
    } else {
        coverage = 0.0;
    }

    const double p_absorb =
        clamp01(coverage * calculate_absorption_prob(e_perp, b_thick));

    return rng.uniform01() < p_absorb;
}

void SurfaceModel::reflect(State& s,
                           const Vec3& norm,
                           const Vec3& tang,
                           RandomEngine& rng) {
    const double p_total = std::sqrt(s.px * s.px + s.py * s.py + s.pz * s.pz);

    const double theta = std::asin(std::sqrt(rng.uniform01()));
    const double phi = 2.0 * constants::kPi * rng.uniform01();

    const double p_n = std::cos(theta);
    const double p_t1 = std::sin(theta) * std::cos(phi);
    const double p_t2 = std::sin(theta) * std::sin(phi);

    const std::vector<double> tang_prime = cross3(norm, tang);

    const double nx = p_n * norm[0] + p_t1 * tang[0] + p_t2 * tang_prime[0];
    const double ny = p_n * norm[1] + p_t1 * tang[1] + p_t2 * tang_prime[1];
    const double nz = p_n * norm[2] + p_t1 * tang[2] + p_t2 * tang_prime[2];

    const double mag = std::sqrt(nx * nx + ny * ny + nz * nz);

    s.px = (nx / mag) * p_total;
    s.py = (ny / mag) * p_total;
    s.pz = (nz / mag) * p_total;
}

} // namespace ucntrap