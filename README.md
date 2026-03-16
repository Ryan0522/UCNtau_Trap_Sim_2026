# UCNtau Trap Sim

A small C++ simulation project for ultra-cold neutron trajectory integration in a magnetic trap.

## Current status

This repository currently includes:

- a planar Halbach field model (4th-order truncated series)
- a simple momentum-based integrator
- a minimal random initial-state generator
- a command-line executable for basic trajectory tests

This is an early working version for development and validation.

## Build

### With CMake

```bash
cmake -S . -B build
cmake --build build --config Debug