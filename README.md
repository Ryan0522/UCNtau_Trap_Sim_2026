# UCNtrap Simulation (v0.1.0)

[cite_start]A high-performance C++ simulation suite designed for Ultra-Cold Neutron (UCN) trajectory integration within the UCNtau magnetic trap [cite: 25-40, 383-392]. [cite_start]This tool supports large-scale ensemble simulations using MPI and provides Python bindings for advanced statistical validation [cite: 115-120, 135-138, 448-452].

## Key Features

* **Physics Models**:
    * [cite_start]**Trap Halbach Field**: Full 3D magnetic field model supporting real-time trap heating/vibrations via coordinate traces [cite: 383-447].
    * [cite_start]**Planar Halbach Field**: 4th-order truncated series model for rapid validation [cite: 76-78, 327-345].
    * [cite_start]**Surface Model**: Quantum mechanical absorption and diffuse reflection logic [cite: 79-81, 346-382].
* **Experimental Logic**:
    * [cite_start]**Three-Phase Tracking**: Dedicated cycles for Cleaning, Holding, and Detection phases to match UCNtau experimental runs [cite: 201-257].
    * [cite_start]**Mechanical Dagger**: Simulated mechanical motion with trapezoidal velocity profiles and acceleration limits [cite: 67-72, 289-326].
* **Performance & Scaling**:
    * [cite_start]**MPI Parallelism**: Efficiently distribute neutron ensembles across multiple CPU ranks with automated result merging [cite: 108-123, 135-138].
    * [cite_start]**Symplectic Integrator**: Energy-conserving 4th-order momentum-based integration for long-term stability [cite: 30-32, 280-288].
    * **Minimalist Progress Bar**: Real-time visual feedback on Rank 0, updating every 5% to minimize I/O overhead.
* **Validation Suite**:
    * [cite_start]**Liouville Test**: Ensures phase-space volume conservation [cite: 477-506].
    * [cite_start]**Distribution Comparison**: Python/SciPy scripts to perform Kolmogorov-Smirnov tests against legacy codebases [cite: 102-103].

## Directory Structure

The build system is configured to output binaries and libraries directly to the project root for streamlined development:

* `./bin/`: Main simulation executables (`ucntrap`).
* `./bin/tests/`: Validation tools (`liouville_test`, `compare_dist_test`).
* `./lib/`: Core physics engine library (`ucntrap_core`).
* [cite_start]`./python/`: Compiled Python module (`ucntrap_py`) and analysis scripts [cite: 102-107, 448-452].

## Build Instructions

### Prerequisites
* CMake 3.20+
* C++17 compliant compiler
* MPI implementation (e.g., OpenMPI or MS-MPI)

### Compilation

To ensure all binaries are placed in the root directories regardless of build configuration:

```bash
# Configure the project
cmake -S . -B build

# Build in Release mode for production performance
cmake --build build --config Release