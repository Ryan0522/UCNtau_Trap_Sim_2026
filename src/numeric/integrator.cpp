#include "ucntrap/numerics/integrator.hpp"
#include "ucntrap/constants.hpp"

#include <array>

namespace ucntrap {

namespace Symp = constants::numeric;

class SymplecticMomentumIntegrator final : public Integrator {
public:
    void step(State& s, double t, double dt, const FieldModel& field) const override {
        constexpr double inv_mass = 1.0 / constants::kMassN;
        
        double local_t = t;

        for (std::size_t n = 0; n < 4; ++n) {
            const Force f = field.force(s, local_t);
            const double kick  = Symp::kSymplecticB[n] * dt;
            const double drift = Symp::kSymplecticA[n] * dt * inv_mass;

            s.px += kick * f.fx;
            s.py += kick * f.fy;
            s.pz += kick * f.fz;

            s.x += drift * s.px;
            s.y += drift * s.py;
            s.z += drift * s.pz;
            
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