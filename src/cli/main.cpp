#include "ucntrap/physics/planar_halbach_field.hpp"
#include "ucntrap/physics/trap_halbach_field.hpp"
#include "ucntrap/numerics/integrator.hpp"
#include "ucntrap/source/random_source.hpp"
#include "ucntrap/config.hpp"
#include "ucntrap/constants.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

namespace {
double kinetic_energy(const ucntrap::State& s)
{
    const double p2 = s.px * s.px + s.py * s.py + s.pz * s.pz;
    return p2 / (2.0 * ucntrap::constants::mass_n);
}
}

int main()
{
    using namespace ucntrap;

    SimulationConfig config;
    config.ntraj = 1;
    config.seed = 123;
    config.dt = 1.0e-3;
    config.heat_mult = 0.0;

    RandomSource source(config);
    State s = source.next();

    s.x = 0.044;
    s.y = -0.322;
    s.z = -1.138;
    s.px = -0.733 * constants::mass_n;
    s.py =  0.118 * constants::mass_n;
    s.pz =  0.003 * constants::mass_n;

    TrapHalbachField field(
        1.35,
        config.heat_mult,
        {}, {}, {}
    );

    const Integrator& integrator = default_integrator();

    double t = 0.0;
    const int steps = 1000;

    std::cout << std::setprecision(12);

    const double U0 = field.potential(s, t);
    const double K0 = kinetic_energy(s);
    const double E0 = K0 + U0;

    auto geom_diag = [](const State& s_now) {
        struct Geom {
            double r_zeta;
            double r;
            bool inside;
        };

        const double x = s_now.x;
        const double y = s_now.y;
        const double z = s_now.z;

        constexpr double KAPPA = 1000.0;
        const double expkx = std::exp(-KAPPA * x);
        const double inv = 1.0 / (1.0 + expkx);
        const double R = 0.5 + 0.5 * inv;
        const double r = 1.0 - 0.5 * inv;

        const double rho = std::sqrt(y * y + z * z);
        const double r_zeta = std::sqrt((rho - R) * (rho - R) + x * x);
        const bool inside = (z < -1.0 && r_zeta < r);

        return Geom{r_zeta, r, inside};
    };

    auto print_brief = [&](const State& s_now, double t_now, const char* tag) {
        const Force f = field.force(s_now, t_now);
        const double U = field.potential(s_now, t_now);
        const double K = kinetic_energy(s_now);
        const double E = K + U;
        const auto g = geom_diag(s_now);

        std::cout << tag
                  << " t=" << t_now
                  << " x=" << s_now.x
                  << " y=" << s_now.y
                  << " z=" << s_now.z
                  << " pz=" << s_now.pz
                  << " fz=" << f.fz
                  << " U=" << U
                  << " K=" << K
                  << " E=" << E
                  << " r_zeta=" << g.r_zeta
                  << " r=" << g.r
                  << " inside=" << g.inside
                  << "\n";
    };

    auto print_full = [&](const State& s_now, double t_now, const char* tag) {
        const Force f = field.force(s_now, t_now);
        const double U = field.potential(s_now, t_now);
        const double K = kinetic_energy(s_now);
        const double E = K + U;
        const auto g = geom_diag(s_now);

        std::cout << tag
                  << " t=" << t_now
                  << " x=" << s_now.x
                  << " y=" << s_now.y
                  << " z=" << s_now.z
                  << " px=" << s_now.px
                  << " py=" << s_now.py
                  << " pz=" << s_now.pz
                  << " fx=" << f.fx
                  << " fy=" << f.fy
                  << " fz=" << f.fz
                  << " U=" << U
                  << " K=" << K
                  << " E=" << E
                  << " r_zeta=" << g.r_zeta
                  << " r=" << g.r
                  << " inside=" << g.inside
                  << "\n";
    };

    print_full(s, t, "initial");

    bool was_inside = geom_diag(s).inside;

    for (int i = 0; i < steps; ++i) {
        integrator.step(s, t, config.dt, field);
        t += config.dt;

        const auto g = geom_diag(s);

        // 每 50 步印一次簡短版
        if (i % 50 == 0) {
            print_brief(s, t, "step");
        }

        // // 快接近邊界時，印完整資訊
        // if (g.inside && (g.r - g.r_zeta) < 0.05) {
        //     print_full(s, t, "near_boundary");
        // }

        // 剛出界的瞬間，印完整資訊
        if (was_inside && !g.inside) {
            print_full(s, t, "exited");
        }

        was_inside = g.inside;
    }

    const double Uf = field.potential(s, t);
    const double Kf = kinetic_energy(s);
    const double Ef = Kf + Uf;

    std::cout << "final:   "
              << s.x << " " << s.y << " " << s.z << " | "
              << s.px << " " << s.py << " " << s.pz << "\n";

    std::cout << "U0 = " << U0 << "\n";
    std::cout << "K0 = " << K0 << "\n";
    std::cout << "E0 = " << E0 << "\n";

    std::cout << "Uf = " << Uf << "\n";
    std::cout << "Kf = " << Kf << "\n";
    std::cout << "Ef = " << Ef << "\n";
    std::cout << "dE = " << (Ef - E0) << "\n";

    return 0;
}