#include "ucntrap/physics/surface_model.hpp"
#include "ucntrap/constants.hpp"

#include <random>


namespace ucntrap {

namespace cs = constants::surface;

namespace {
    using ComplexMat = std::array<std::complex<double>, 4>;
    
    ComplexMat matmul(const ComplexMat& a, const ComplexMat& b) {
        return {
            a[0]*b[0] + a[1]*b[2], a[0]*b[1] + a[1]*b[3],
            a[2]*b[0] + a[3]*b[2], a[2]*b[1] + a[3]*b[3]
        };
    }

    std::complex<double> k_wave(double e_perp, std::complex<double> u) {
        return std::sqrt((2.0 * constants::kMassN / (constants::kHbar * constants::kHbar)) * (e_perp - u));
    }

    ComplexMat m_matrix(std::complex<double> kn, std::complex<double> kn_minus_1, double z) {
        std::complex<double> gamma = kn_minus_1 / kn;
        std::complex<double> i(0.0, 1.0);
        
        ComplexMat res;
        res[0] = 0.5 * (1.0 + gamma) * std::exp(i * (kn_minus_1 - kn) * z);
        res[1] = 0.5 * (1.0 - gamma) * std::exp(-i * (kn_minus_1 + kn) * z);
        res[2] = 0.5 * (1.0 - gamma) * std::exp(i * (kn_minus_1 + kn) * z);
        res[3] = 0.5 * (1.0 + gamma) * std::exp(-i * (kn_minus_1 - kn) * z);
        return res;
    }
}

std::vector<std::complex<double>> matmul(const std::vector<std::complex<double>>& a, 
                                        const std::vector<std::complex<double>>& b) {
    return { a[0]*b[0] + a[1]*b[2], a[0]*b[1] + a[1]*b[3],
             a[2]*b[0] + a[3]*b[2], a[2]*b[1] + a[3]*b[3] };
}

double SurfaceModel::calculate_absorption_prob(double e_perp, double b_thick) {
    // 1. Calculate potential
    double v_boron = (2.0 * constants::kPi * (constants::kHbar * constants::kHbar) / constants::kMassN) * cs::kABoron * cs::kNBoron;
    double w_boron = (constants::kHbar / 2.0) * cs::kNBoron * cs::kSigmaBoron;
    
    const double v_zns = (2.0 * constants::kPi * (constants::kHbar * constants::kHbar) / constants::kMassN) * cs::kAZinc * cs::kNZinc 
                         + (2.0 * constants::kPi * (constants::kHbar * constants::kHbar) / constants::kMassN) * cs::kASulfur * cs::kNSulfur;
    const double w_zns = (constants::kHbar / 2.0) * (cs::kNZinc * cs::kSigmaZinc + cs::kNSulfur * cs::kSigmaSulfur);

    std::vector<std::complex<double>> potentials = {
        std::complex<double>(0.0, 0.0),             // vacuum
        std::complex<double>(v_boron, -w_boron),    // Boron
        std::complex<double>(v_zns, -w_zns)         // Zns
    };

    std::vector<double> z_boundaries = { 0.0, b_thick_nm * 1e-9, 10000e-9 }; // in meters

    ComplexMat m_bar = { 1.0, 0.0, 0.0, 1.0 }; // Identity matrix
    for (int i = potentials.size() - 1; i > 0; --i) {
        auto kn = k_wave(e_perp, potentials[i]);
        auto kn_m1 = k_wave(e_perp, potentials[i-1]);
        m_bar = matmul(m_bar, m_matrix(kn, kn_m1, z_boundaries[i-1]));
    }
    
    // Calculate Reflection R = -M21 / M22，Absorption P = 1 - |R|^2
    std::complex<double> r = -m_bar[2] / m_bar[3];
    return 1.0 - std::norm(r);
}

bool SurfaceModel::check_absorption(double e_perp, double b_thick_nm, double x, double z, double z_off) {
    // Calculate boron layer effective coverage based on zeta
    double z_rel = std::abs(z - z_off);
    double zeta = (x > 0.0) ? (0.5 - std::sqrt(x * x + std::pow(z_rel - 1.0, 2))) 
                            : (1.0 - std::sqrt(x * x + std::pow(z_rel - 0.5, 2)));

    double coverage = (zeta > constants::kZetaCut) ? 1.0 : (zeta / constants::kZetaCut);
    
    static std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng) < (coverage * calculate_absorption_prob(e_perp, b_thick_nm));
}

void SurfaceModel::reflect(State& s, const std::vector<double>& norm, const std::vector<double>& tang) {
    static std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    double p_total = std::sqrt(s.px*s.px + s.py*s.py + s.pz*s.pz);
    double theta = std::asin(std::sqrt(dist(rng))); 
    double phi = 2.0 * constants::kPi * dist(rng);

    double p_n = std::cos(theta);
    double p_t1 = std::sin(theta) * std::cos(phi);
    double p_t2 = std::sin(theta) * std::sin(phi);

    std::vector<double> tang_prime = {
        norm[1] * tang[2] - norm[2] * tang[1],
        norm[2] * tang[0] - norm[0] * tang[2],
        norm[0] * tang[1] - norm[1] * tang[0]
    };

    double nx = p_n * norm[0] + p_t1 * tang[0] + p_t2 * tang_prime[0];
    double ny = p_n * norm[1] + p_t1 * tang[1] + p_t2 * tang_prime[1];
    double nz = p_n * norm[2] + p_t1 * tang[2] + p_t2 * tang_prime[2];

    double mag = std::sqrt(nx*nx + ny*ny + nz*nz);
    s.px = (nx / mag) * p_total;
    s.py = (ny / mag) * p_total;
    s.pz = (nz / mag) * p_total;
}

} // namespace ucntrap