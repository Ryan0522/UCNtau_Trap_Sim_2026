#include "ucntrap/runner.hpp"
#include "ucntrap/config.hpp"

#include <mpi.h>
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    try {
        // read YAML
        ucntrap::SimulationConfig config = ucntrap::load_config("config/config.yaml");

        // simulate
        ucntrap::Runner runner(config);
        runner.run();

    } catch (const std::exception& e) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            std::cerr << "Fatal Error: " << e.what() << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}