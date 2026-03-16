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
    config.dt = 1.0e-5;
    config.heat_mult = 0.0;

    RandomSource source(config);
    State s = source.next();

    s.x = -0.0075;
    s.y = 0.0022;
    s.z = -1.10;
    s.px = 7.9e-28;
    s.py = -1.03e-27;
    s.pz = -1.02e-27;

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

    std::cout << "initial: "
              << s.x << " " << s.y << " " << s.z << " | "
              << s.px << " " << s.py << " " << s.pz << "\n";
    std::cout << "U0 = " << U0 << "\n";
    std::cout << "K0 = " << K0 << "\n";
    std::cout << "E0 = " << E0 << "\n";

    for (int i = 0; i < steps; ++i) {
        integrator.step(s, t, config.dt, field);
        t += config.dt;

        if (i % 100 == 0) {
            const double U = field.potential(s, t);
            const double K = kinetic_energy(s);
            const double E = K + U;

            std::cout << "step " << i
                      << "  x=" << s.x
                      << "  y=" << s.y
                      << "  z=" << s.z
                      << "  E=" << E
                      << "  dE=" << (E - E0)
                      << "\n";
        }
    }

    const double Uf = field.potential(s, t);
    const double Kf = kinetic_energy(s);
    const double Ef = Kf + Uf;

    std::cout << "final:   "
              << s.x << " " << s.y << " " << s.z << " | "
              << s.px << " " << s.py << " " << s.pz << "\n";

    std::cout << "Uf = " << Uf << "\n";
    std::cout << "Kf = " << Kf << "\n";
    std::cout << "Ef = " << Ef << "\n";
    std::cout << "dE = " << (Ef - E0) << "\n";

    return 0;
}