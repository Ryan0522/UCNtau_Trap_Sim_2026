#include "ucntrap/runner.hpp"
#include "ucntrap/source/pentrack_reader.hpp"
#include "ucntrap/config.hpp"

#include <mpi.h>
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    try {
        ucntrap::SimulationConfig config;
        config.dt = 0.001;
        config.ntraj = 1000;
        config.source_type = "pentrack";
        config.neutron_init_file = "./data/neutrons_init_test.out";
        config.field_model = "trap";
        config.output_prefix = "./tests/files/modern_dist";
        config.x_trace_file = "./data/xvals.bin";
        config.y_trace_file = "./data/yvals.bin";
        config.z_trace_file = "./data/zvals.bin";
        config.hold_time = 20.0;
        config.heat_mult = 0.0;

        ucntrap::Runner runner(config);
        runner.run();

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            std::cout << "Distribution simulation complete. Results saved with prefix: "
                        << config.output_prefix << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error during distribution test: " << e.what() << std::endl;
    }

    MPI_Finalize();
    return 0;
}