#include "ucntrap/physics/planar_halbach_field.hpp"
#include "ucntrap/constants.hpp"

#include <cmath>
#include <algorithm>

namespace ucntrap {

namespace {
    
struct FieldEval {
    double bx = 0.0;
    double bz = 0.0;
    double bmag = 0.0;
    double dBxdx = 0.0;
    double dBzdx = 0.0;
    double dBxdz = 0.0;
    double dBzdz = 0.0;
};

FieldEval eval_planar_halbach(double x, double z,
                                double b_rem,
                                double mag_space,
                                double mag_thick,
                                int n_terms)
{
    FieldEval out{};

    const double A = 4.0 * b_rem / (constants::kPi * std::sqrt(2.0));
    const double z_shift = z + 1.5;   // 你原本的 offset
    const double zz = z_shift;

    for (int n = 1; n <= n_terms; ++n) {
        const double odd = 4.0 * n - 3.0;
        const double k = 2.0 * constants::kPi * odd / mag_space;
        const double sign = (n % 2 == 0) ? 1.0 : -1.0;
        const double amp = sign / odd
                            * (1.0 - std::exp(-k * mag_thick))
                            * std::exp(-k * zz);

        const double c = std::cos(k * x);
        const double s = std::sin(k * x);

        const double bx_term = amp * s;
        const double bz_term = amp * c;

        out.bx += bx_term;
        out.bz += bz_term;

        out.dBxdx += k * bz_term;
        out.dBzdx -= k * bx_term;
        out.dBxdz -= k * bx_term;
        out.dBzdz -= k * bz_term;
    }

    out.bx *= A;
    out.bz *= A;
    out.dBxdx *= A;
    out.dBzdx *= A;
    out.dBxdz *= A;
    out.dBzdz *= A;

    out.bmag = std::sqrt(out.bx * out.bx + out.bz * out.bz);
    return out;
}

} // namespace

PlanarHalbachField::PlanarHalbachField(double b_rem, double mag_space,
                                       double mag_thick, int n_terms)
    : b_rem_(b_rem), mag_space_(mag_space), mag_thick_(mag_thick), n_terms_(n_terms) {}

Force PlanarHalbachField::force(const State& s, double /*t*/) const
{
    const auto f = eval_planar_halbach(s.x, s.z, b_rem_, mag_space_, mag_thick_, n_terms_);
    const double bmag = std::max(f.bmag, constants::kEpsilon);

    Force out{};
    out.fx = -constants::kMuN * (f.bx * f.dBxdx + f.bz * f.dBzdx) / bmag;
    out.fy = 0.0;
    out.fz = -constants::kMuN * (f.bx * f.dBxdz + f.bz * f.dBzdz) / bmag
           - constants::kMassN * constants::kEarthG;
    return out;
}

double PlanarHalbachField::potential(const State& s, double /*t*/) const
{
    const auto f = eval_planar_halbach(s.x, s.z, b_rem_, mag_space_, mag_thick_, n_terms_);
    const double z_grav = s.z + 1.5;
    return constants::kMuN * f.bmag + constants::kMassN * constants::kEarthG * z_grav;
}

} // namespace ucntrap