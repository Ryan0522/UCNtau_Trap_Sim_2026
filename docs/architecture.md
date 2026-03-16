

MPI Intialize

Config

readTrace

RNG Intialize

Particles Initialize
    Energy
    Potential
        Halbach 4th Order Expansion
    Angles

Print Flag

Integration
    Potential
    Decay Time
    Starting Time (Settling + Clean) 
        # Found issue: should not clean when settling
    Steps (Check clean and after)
        Integration Step
            Find force
                Define local coordinates
                Find F_g
                Find F_B
            4th roder
        Check decay
        Check clean
        Check energy
        Check crossing dagger plane (Hit -> Absorb? -> Reflect)
            Check dagger hit
            Check house hit low
            Check house hit high

Write to File