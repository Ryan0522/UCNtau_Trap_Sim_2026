#include "ucntrap/numerics/integrator.hpp"
#include "ucntrap/constants.hpp"

#include <array>

namespace ucntrap {

namespace Symp = constants::numeric;

class SymplecticMomentumIntegrator final : public Integrator {
public:
    void step(State& s, double t, double dt, const FieldModel& field) const override {
        double local_t = t;

        for (std::size_t n = 0; n < 4; ++n) {
            const Force f = field.force(s, local_t);

            // kick
            s.px += Symp::kSymplecticB[n] * f.fx * dt;
            s.py += Symp::kSymplecticB[n] * f.fy * dt;
            s.pz += Symp::kSymplecticB[n] * f.fz * dt;

            // drift
            s.x += Symp::kSymplecticA[n] * s.px * dt / constants::kMassN;
            s.y += Symp::kSymplecticA[n] * s.py * dt / constants::kMassN;
            s.z += Symp::kSymplecticA[n] * s.pz * dt / constants::kMassN;

            local_t += Symp::kSymplecticA[n] * dt;
        }
    }
};

const Integrator& default_integrator()
{
    static SymplecticMomentumIntegrator inst;
    return inst;
}

} // namespace ucntrap