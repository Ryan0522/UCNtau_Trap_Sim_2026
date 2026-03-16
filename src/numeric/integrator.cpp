#include "ucntrap/numerics/integrator.hpp"
#include "ucntrap/constants.hpp"

#include <array>

namespace ucntrap {
namespace {
inline constexpr std::array<double, 4> kA = {
    0.5153528374311229364,
   -0.0857820194129736460,
    0.4415830236164665242,
    0.1288461583653841854
};

inline constexpr std::array<double, 4> kB = {
    0.1344961992774310892,
   -0.2248198030794208058,
    0.7563200005156682911,
    0.3340036032863214255
};
} // namespace

class SymplecticMomentumIntegrator final : public Integrator {
public:
    void step(State& s, double t, double dt, const FieldModel& field) const override {
        double local_t = t;

        for (std::size_t n = 0; n < 4; ++n) {
            const Force f = field.force(s, local_t);

            // kick
            s.px += kB[n] * f.fx * dt;
            s.py += kB[n] * f.fy * dt;
            s.pz += kB[n] * f.fz * dt;

            // drift
            s.x += kA[n] * s.px * dt / constants::mass_n;
            s.y += kA[n] * s.py * dt / constants::mass_n;
            s.z += kA[n] * s.pz * dt / constants::mass_n;

            local_t += kA[n] * dt;
        }
    }
};

const Integrator& default_integrator()
{
    static SymplecticMomentumIntegrator inst;
    return inst;
}

} // namespace ucntrap