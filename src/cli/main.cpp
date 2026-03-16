#include "ucntrap/runner.hpp"
#include "ucntrap/config.hpp"

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void merge_csv_files(const std::string& prefix, int num_ranks) {
    std::ofstream final_out(prefix + ".csv");
    bool header_copied = false;

    for (int i = 0; i < num_ranks; ++i) {
        std::string rank_file = prefix + "_rank" + std::to_string(i) + ".csv";
        std::ifstream input(rank_file);

        if (!input.is_open()) continue;

        std::string line;
        bool is_first_line = true;
        while (std::getline(input, line)) {
            if (is_first_line) {
                if (!header_copied) {
                    final_out << line << "\n";
                    header_copied = true;
                }
                is_first_line = false;
                continue;
            }
            final_out << line << "\n";
        }
        input.close();
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    try {
        // read YAML
        ucntrap::SimulationConfig config = ucntrap::load_config("config/config.yaml");

        // simulate
        ucntrap::Runner runner(config);
        runner.run();

        // Make sure all MPI are done
        MPI_Barrier(MPI_COMM_WORLD);

        // Using rank = 0 to merge csvs
        if (rank == 0) {
            std::cout << "Merging result files..." << std::endl;
            merge_csv_files(config.output_prefix, size);
            std::cout << "Merge complete: " << config.output_prefix << ".csv" << std::endl;
        }

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