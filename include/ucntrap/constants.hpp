#pragma once
#include <array>

namespace ucntrap::constants {

// --- 1. basic math and physics constants ---
inline constexpr double kPi      = 3.14159265358979323846;
inline constexpr double kEpsilon = 1.0e-9;
inline constexpr double kMuN     = 9.6623647e-27;
inline constexpr double kMassN   = 1.674927471e-27;
inline constexpr double kEarthG  = 9.80665e0;
inline constexpr double kHbar    = 1.054571800e-34;
inline constexpr double kJToEv   = 6.2415091e27;
inline constexpr double kMinU    = -2.390352484438862e-26;

// --- 2. axes scaling (for Lyapunov and/or momentum) ---
inline constexpr double kXScale = 0.2;
inline constexpr double kYScale = 0.35;
inline constexpr double kZScale = 1.34;
inline constexpr double kPScale = 1.25e-27;

// --- 3. UCNtau Apparatus ---
namespace trap {
    // Halbach Array Parameters
    inline constexpr double kMagSpace  = 0.05114;
    inline constexpr double kMagThick  = 0.0254;
    inline constexpr double kKappa     = 1000.0;
    inline constexpr double kBHold     = 0.005;
    inline constexpr int    kNSumTerms = 3;

    // Oscillation sampling (For "Shift" function using Trace)
    inline constexpr double kSampleDt = 0.0004;
}

// --- 4. Symplectic Integrator Parameters ---
namespace numeric {
    inline constexpr std::array<double, 4> kSymplecticA = {
        0.5153528374311229364,
       -0.0857820194129736460,
        0.4415830236164665242,
        0.1288461583653841854
    };
    inline constexpr std::array<double, 4> kSymplecticB = {
        0.1344961992774310892,
       -0.2248198030794208058,
        0.7563200005156682911,
        0.3340036032863214255
    };
}

} // namespace ucntrap::constants