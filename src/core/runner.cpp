#include "ucntrap/runner.hpp"
#include "ucntrap/config.hpp"

#include <mpi.h>
#include <iostream>
#include <memory>

namespace ucntrap {

Runner::Runner(SimulationConfig config) : config_(config) {}

int Runner::run() const {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 1. Calculate neutron counts in this rank
    size_t total_ntraj = config_.ntraj;
    size_t my_ntraj = total_ntraj / size;
    if (rank == size - 1) {
        my_ntraj += (total_ntraj % size);
    }

    if (rank == 0) {
        std::cout << "Starting Simulation with " << size << " ranks." << std::endl;
        std::cout << "Total Trajectories: " << total_ntraj << std::endl;
    }

    return 0;
}

} // namespace ucntrap