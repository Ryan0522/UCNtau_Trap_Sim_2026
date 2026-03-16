#include "ucntrap/numerics/integrator.hpp"
#include "ucntrap/constants.hpp"

namespace ucntrap {

class SymplecticMomentumIntegrator final : public Integrator {
public:
    void step(State& s, double t, double dt, const FieldModel& field) const override {
        // kick: p_{n+1/2} = p_n + (dt/2) F(x_n)
        const Force f0 = field.force(s, t);

        s.px += 0.5 * dt * f0.fx;
        s.py += 0.5 * dt * f0.fy;
        s.pz += 0.5 * dt * f0.fz;

        // drift: x_{n+1} = x_n + dt * p_{n+1/2} / m
        s.x += dt * s.px / constants::mass_n;
        s.y += dt * s.py / constants::mass_n;
        s.z += dt * s.pz / constants::mass_n;

        // kick: p_{n+1} = p_{n+1/2} + (dt/2) F(x_{n+1})
        const Force f1 = field.force(s, t + dt);

        s.px += 0.5 * dt * f1.fx;
        s.py += 0.5 * dt * f1.fy;
        s.pz += 0.5 * dt * f1.fz;
    }
};

const Integrator& default_integrator()
{
    static SymplecticMomentumIntegrator inst;
    return inst;
}

} // namespace ucntrap