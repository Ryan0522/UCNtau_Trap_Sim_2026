#include "ucntrap/physics/trap_halbach_field.hpp"
#include "ucntrap/constants.hpp"

#include <cmath>

namespace ucntrap {

namespace ct = constants::trap;

TrapHalbachField::TrapHalbachField(double heat_mult,
                     const std::vector<double>& tx,
                     const std::vector<double>& ty,
                     const std::vector<double>& tz)
    : heat_mult_(heat_mult), tx_(tx), ty_(ty), tz_(tz)
{
    A_ = 4.0 * ct::kBRem / (constants::kPi * std::sqrt(2.0));

    for (int n = 1; n <= ct::kNSumTerms; ++n) {
        const int idx = n - 1;
        const double odd = 4.0 * n - 3.0;
        const double sign = (n % 2 == 0) ? 1.0 : -1.0;

        k_n_[idx] = 2.0 * constants::kPi * odd / ct::kMagSpace;
        amp_prefactor_[idx] = sign / odd * (1.0 - std::exp(-k_n_[idx] * ct::kMagThick));
    }
}

void TrapHalbachField::get_shifted_coords(const State& s, double t, double& x, double& y, double& z) const {
    
    if (heat_mult_ == 0.0 || tx_.empty() || ty_.empty() || tz_.empty()) {
        x = s.x;
        y = s.y;
        z = s.z;
        return;
    }
    
    int num = tx_.size();
    int iLowRaw = static_cast<int>(t / ct::kSampleDt);
    double frac = (t - iLowRaw * ct::kSampleDt) / ct::kSampleDt;

    int iLow = (iLowRaw % num + num) % num;
    int iHi  = ((iLowRaw + 1) % num + num) % num;

    x = s.x + heat_mult_ * (tx_[iLow] + frac * (tx_[iHi] - tx_[iLow]));
    y = s.y + heat_mult_ * (ty_[iLow] + frac * (ty_[iHi] - ty_[iLow]));
    z = s.z + heat_mult_ * (tz_[iLow] + frac * (tz_[iHi] - tz_[iLow]));
}

Force TrapHalbachField::force(const State& s, double t) const {
    double x, y, z;
    get_shifted_coords(s, t, x, y, z);

    const double expkx = std::exp(-ct::kKappa * x);
    const double inv = 1.0 / (1.0 + expkx);

    const double R = 0.5 + 0.5 * inv;
    const double r = 1.0 - 0.5 * inv;
    const double Rprime = 0.5 * ct::kKappa * (1.0 - inv) * inv;
    const double rprime = -Rprime;

    const double rho = std::sqrt(y * y + z * z);
    const double safe_rho = std::max(rho, constants::kEpsilon);
    const double delta = safe_rho - R;
    const double root = std::sqrt(x * x + delta * delta);
    const double r_zeta = root;

    Force f_out{NAN, NAN, NAN};

    if (z < -1.0 && r_zeta < r) {
        const double eta = r * std::atan(x / (safe_rho - R));
        const double zeta = r - root;

        double sum_cos = 0.0;
        double sum_sin = 0.0;
        double sum_k_cos = 0.0;
        double sum_k_sin = 0.0;

        for (int idx = 0; idx < ct::kNSumTerms; ++idx) {
            const double k = k_n_[idx];
            const double amp = amp_prefactor_[idx] * std::exp(-k * zeta);

            const double c = std::cos(k * eta);
            const double s_eta = std::sin(k * eta);

            const double cos_term = amp * c;
            const double sin_term = amp * s_eta;

            sum_cos += cos_term;
            sum_sin += sin_term;
            sum_k_cos += k * cos_term;
            sum_k_sin += k * sin_term;
        }

        const double b_zeta = A_ * sum_cos;
        const double b_eta  = A_ * sum_sin;
        const double b_hold = ct::kBHold * (r + R) / safe_rho;

        const double b_tot = std::sqrt(
            b_zeta * b_zeta +
            b_eta  * b_eta +
            b_hold * b_hold
        );
        const double safe_b_tot = std::max(b_tot, constants::kEpsilon);

        // derivatives wrt (zeta, eta)
        const double d_BZeta_0 = -A_ * sum_k_cos;
        const double d_BZeta_1 = -A_ * sum_k_sin;

        const double d_BEta_0 = -A_ * sum_k_sin;
        const double d_BEta_1 =  A_ * sum_k_cos;

        // derivatives wrt (x, y, z)
        const double d_Bh_x = 0.0;
        const double d_Bh_y = -ct::kBHold * (r + R) * y / (safe_rho * safe_rho * safe_rho);
        const double d_Bh_z = -ct::kBHold * (r + R) * z / (safe_rho * safe_rho * safe_rho);

        const double denom_root = std::max(root, constants::kEpsilon);
        const double R_minus_rho = R - safe_rho;
        const double denom_eta_1 = x * x + y * y + z * z + R * R - 2.0 * R * safe_rho;
        const double safe_denom_eta_1 = std::max(std::abs(denom_eta_1), constants::kEpsilon)
                                      * (denom_eta_1 < 0.0 ? -1.0 : 1.0);
        const double denom_eta_2 = (R_minus_rho * R_minus_rho)
                                 * (1.0 + x * x / std::max(R_minus_rho * R_minus_rho, constants::kEpsilon));

        const double d_Zeta_x =
            rprime + (Rprime * (safe_rho - R) - x) / denom_root;
        const double d_Zeta_y =
            -y * (safe_rho - R) / (safe_rho * denom_root);
        const double d_Zeta_z =
            -z * (safe_rho - R) / (safe_rho * denom_root);

        const double d_Eta_x =
            -rprime * std::atan(x / (R - safe_rho))
            - (r * (R - safe_rho - x * Rprime)) / safe_denom_eta_1;

        const double d_Eta_y =
            -r * x * y / (safe_rho * std::max(denom_eta_2, constants::kEpsilon));

        const double d_Eta_z =
            -r * x * z / (safe_rho * std::max(denom_eta_2, constants::kEpsilon));

        f_out.fx = -constants::kMuN * (1.0 / safe_b_tot) * (
            b_zeta * (d_BZeta_0 * d_Zeta_x + d_BZeta_1 * d_Eta_x) +
            b_eta  * (d_BEta_0  * d_Zeta_x + d_BEta_1  * d_Eta_x) +
            b_hold * d_Bh_x
        );

        f_out.fy = -constants::kMuN * (1.0 / safe_b_tot) * (
            b_zeta * (d_BZeta_0 * d_Zeta_y + d_BZeta_1 * d_Eta_y) +
            b_eta  * (d_BEta_0  * d_Zeta_y + d_BEta_1  * d_Eta_y) +
            b_hold * d_Bh_y
        );

        f_out.fz = -constants::kMuN * (1.0 / safe_b_tot) * (
            b_zeta * (d_BZeta_0 * d_Zeta_z + d_BZeta_1 * d_Eta_z) +
            b_eta  * (d_BEta_0  * d_Zeta_z + d_BEta_1  * d_Eta_z) +
            b_hold * d_Bh_z
        );

        f_out.fz -= constants::kMassN * constants::kEarthG;
    }

    return f_out;
}

double TrapHalbachField::potential(const State& s, double t) const {
    double x, y, z;
    get_shifted_coords(s, t, x, y, z);

    // gravity uses unshifted z, matching legacy code
    const double z_grav = s.z;

    const double expkx = std::exp(-ct::kKappa * x);
    const double inv = 1.0 / (1.0 + expkx);

    const double R = 0.5 + 0.5 * inv;
    const double r = 1.0 - 0.5 * inv;

    const double rho = std::sqrt(y * y + z * z);
    const double safe_rho = std::max(rho, constants::kEpsilon);
    const double r_zeta = std::sqrt((safe_rho - R) * (safe_rho - R) + x * x);

    if (!(z < -1.0 && r_zeta < r)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double eta = r * std::atan(x / (safe_rho - R));
    const double zeta = r - std::sqrt(x * x + (safe_rho - R) * (safe_rho - R));

    double sum_cos = 0.0;
    double sum_sin = 0.0;

    for (int idx = 0; idx < ct::kNSumTerms; ++idx) {
        const double k = k_n_[idx];
        const double amp = amp_prefactor_[idx] * std::exp(-k * zeta);

        sum_cos += amp * std::cos(k * eta);
        sum_sin += amp * std::sin(k * eta);
    }

    const double b_zeta = A_ * sum_cos;
    const double b_eta  = A_ * sum_sin;
    const double b_hold = ct::kBHold * (r + R) / safe_rho;
    const double b_tot = std::sqrt(
        b_zeta * b_zeta +
        b_eta  * b_eta +
        b_hold * b_hold
    );

    return constants::kMuN * b_tot + constants::kMassN * constants::kEarthG * z_grav;
}

} // namespace ucntrap 