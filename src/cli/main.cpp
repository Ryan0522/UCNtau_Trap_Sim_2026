#include "ucntrap/physics/field_model.hpp"
#include "ucntrap/numerics/integrator.hpp"
#include "ucntrap/source/random_source.hpp"

#include <iostream>
#include <iomanip>

int main() 
{
    using namespace ucntrap;

    SimulationConfig config;
    config.ntraj = 1;
    config.seed = 123;
    config.dt = 1.0e-5;

    RandomSource source(config);
    State s = source.next();

    PlanarHalbachField field(
        1.32, // B_rem
        0.05, // MAG_SPACE
        0.025, // MAG_THICK
        4 // N_TERMS
    );

    const Integrator& integrator = default_integrator();

    double t = 0.0;
    const int step = 1000;

    std::cout << std::setprecision(10);
    std::cout << "initial: "
              << s.x << " " << s.y << " " << s.z << " | "
              << s.px << " " << s.py << " " << s.pz << "\n";

    for (int i = 0; i < steps; ++i) {
        integrator.step(s, t, config.dt, field);
        t += config.dt;
    }

    std::cout << "final:   "
              << s.x << " " << s.y << " " << s.z << " | "
              << s.px << " " << s.py << " " << s.pz << "\n";

    std::cout << "U_final = " << field.potential(s, t) << "\n";
    return 0;
}