#pragma once
#include "ucntrap/physics/field_model.hpp"
#include "ucntrap/constants.hpp"
#include "ucntrap/numerics/lut.hpp"
#include <array>
#include <vector>

namespace ucntrap {

class TrapHalbachField final : public FieldModel {
public:
    TrapHalbachField(double heat_mult,
                     const std::vector<double>& tx,
                     const std::vector<double>& ty,
                     const std::vector<double>& tz);
    
    Force force(const State& s, double t) const override;
    double potential(const State& s, double t) const override;

private:
    static constexpr int kTerms = constants::trap::kNSumTerms;
    double heat_mult_;
    const std::vector<double> tx_, ty_, tz_;

    LinearLUT sin_, cos_, exp_;

    double A_ = 0.0;
    std::array<double, kTerms> k_n_{};
    std::array<double, kTerms> amp_prefactor_{};

    void get_shifted_coords(const State& s, double t, double& x, double& y, double& z) const;
};

} // namespace ucntrap