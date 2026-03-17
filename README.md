# UCNtrap Simulation (v0.1.0)

A high-performance C++ simulation suite designed for Ultra-Cold Neutron (UCN) trajectory integration within the UCNtau magnetic trap. This tool supports large-scale ensemble simulations using MPI and provides Python bindings for advanced statistical validation and data analysis.

## Key Features

* **Physics Models**:
    * **Trap Halbach Field**: Full 3D magnetic field model supporting real-time trap heating/vibrations via coordinate traces.
    * **Planar Halbach Field**: 4th-order truncated series model for rapid validation and testing.
    * **Surface Model**: Quantum mechanical absorption and diffuse reflection logic based on material properties.
* **Experimental Logic**:
    * **Three-Phase Tracking**: Dedicated cycles for Cleaning, Holding, and Detection phases to accurately replicate UCNtau experimental runs.
    * **Mechanical Dagger**: Simulated mechanical motion with trapezoidal velocity profiles and configurable acceleration limits.
* **Performance & Scaling**:
    * **MPI Parallelism**: Efficiently distribute neutron ensembles across multiple CPU ranks with automated result merging.
    * **Symplectic Integrator**: Energy-conserving 4th-order momentum-based integration for long-term orbital stability.
* **Validation Suite**:
    * **Trajectory Test**: Built-in tool to ensure trajectories are consistent until Lyapunov time.
    * **Liouville Test**: Built-in tool to ensure phase-space volume conservation.
    * **Distribution Comparison**: Python/SciPy scripts to perform Kolmogorov-Smirnov tests against legacy codebases.

## Directory Structure

The build system is configured to output binaries and libraries directly to the project root for streamlined development and access:

* `./bin/`: Main simulation executables (`ucntrap`).
* `./bin/tests/`: Validation and stability tools (`traj_overlap_test`, `liouville_test`, `compare_dist_test`).
* `./lib/`: Core physics engine static/shared libraries.
* `./python/`: Compiled Python module (`ucntrap_py`) and analysis scripts.

## Build Instructions

### Prerequisites
* CMake 3.20+
* C++17 compliant compiler (GCC, Clang, or MSVC)
* MPI implementation (e.g., OpenMPI, MPICH, or MS-MPI)

### Compilation

To ensure all binaries are placed in the root directories regardless of build configuration (Debug/Release):

```bash
# Configure the project
cmake -S . -B build

# Build in Release mode for production-grade performance
cmake --build build --config Release
```

## Usage

**1. Configuration**

Modify `config/config.yaml` to set experimental timings (cleaning/hold times), dagger motion schedules, and core physics/experiment parameters.

**2. Running Simulations**

Execute the MPI-enabled runner from the root directory:

```bash
mpirun -np 4 ./bin/ucntrap
```

**3. Statistical Analysis**

Post-analysis for simulation result in Python.
