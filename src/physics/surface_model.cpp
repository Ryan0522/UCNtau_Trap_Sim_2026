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
        return std::sqrt((2.0 * un::kMassN / (un::kHbar * un::kHbar)) * (e_perp - u));
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
    double v_boron = (2.0 * un::kPi * (un::kHbar * un::kHbar) / un::kMassN) * un::surface::kABoron * un::surface::kNBoron;
    double w_boron = (un::kHbar / 2.0) * un::surface::kNBoron * un::surface::kSigmaBoron;
    
    const double v_zns = (2.0 * un::kPi * (un::kHbar * un::kHbar) / un::kMassN) * un::surface::kAZinc * un::surface::kNZinc 
                         + (2.0 * un::kPi * (un::kHbar * un::kHbar) / un::kMassN) * un::surface::kASulfur * un::surface::kNSulfur;
    const double w_zns = (un::kHbar / 2.0) * (un::surface::kNZinc * un::surface::kSigmaZinc + un::surface::kNSulfur * un::surface::kSigmaSulfur);

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
    double zeta;
    if (x > 0.0) {
        zeta = 0.5 - std::sqrt(x * x + std::pow(std::abs(z - z_off) - 1.0, 2));
    } else {
        zeta = 1.0 - std::sqrt(x * x + std::pow(std::abs(z - z_off) - 0.5, 2));
    }

    double coverage = (zeta > un::kZetaCut) ? 1.0 : (zeta / un::kZetaCut);
    
    static std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    return dist(rng) < (coverage * calculate_absorption_prob(e_perp, b_thick_nm));
}

void SurfaceModel::reflect(State& s, const std::vector<double>& norm) {
    static std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    double p_total = std::sqrt(s.px*s.px + s.py*s.py + s.pz*s.pz);
    
    double theta = std::asin(std::sqrt(dist(rng))); 
    double phi = 2.0 * un::kPi * dist(rng);

    double p_normal = p_total * std::cos(theta);
    double p_tangent1 = p_total * std::sin(theta) * std::cos(phi);
    double p_tangent2 = p_total * std::sin(theta) * std::sin(phi);

    s.py = (norm[1] > 0) ? std::abs(p_normal) : -std::abs(p_normal);
    s.px = p_tangent1;
    s.pz = p_tangent2;
}

} // namespace ucntrap