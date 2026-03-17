#include "ucntrap/runner.hpp"
#include "ucntrap/experiment/production_tracker.hpp"
#include "ucntrap/physics/trap_halbach_field.hpp"
#include "ucntrap/physics/planar_halbach_field.hpp"
#include "ucntrap/numerics/integrator.hpp"
#include "ucntrap/source/pentrack_reader.hpp"
#include "ucntrap/source/random_source.hpp"
#include "ucntrap/io/trace_loader.hpp"
#include "ucntrap/io/result_writer.hpp"

#include <mpi.h>
#include <iostream>
#include <memory>

namespace ucntrap {

Runner::Runner(SimulationConfig config) : config_(config) {}

int Runner::run() const {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    RandomEngine rng(config_.seed + rank);

    // 1. Calculate neutron counts in this rank
    size_t total_ntraj = config_.ntraj;
    size_t my_ntraj = total_ntraj / size;
    size_t remainder = total_ntraj % size;

    if (rank == size - 1) {
        my_ntraj += remainder;
    }

    // 2. Calculate offset in reading PENTrack source
    size_t global_offset = 1 + (rank * (total_ntraj / size));

    if (rank == 0) {
        print_config(config_, size, total_ntraj);
    }

    // 3. Initialize source
    std::unique_ptr<InitialConditionSource> source;
    if (config_.source_type == "pentrack") {
        source = std::make_unique<PenTrackReader>(config_.neutron_init_file, my_ntraj, global_offset);
    } else {
        source = std::make_unique<RandomSource>(config_);
    }

    // 4. Initialize field
    std::unique_ptr<FieldModel> field;
    if (config_.field_model == "trap") {
        Trace trace_data = load_trace(config_.x_trace_file,
                                        config_.y_trace_file,
                                        config_.z_trace_file);

        field = std::make_unique<TrapHalbachField>(
            config_.heat_mult,
            trace_data.x,
            trace_data.y,
            trace_data.z
        );
    } else {
        field = std::make_unique<PlanarHalbachField>(1.35, 0.05114, 0.0254, 3);
    }

    // 5. Initialize integrators and trackers
    const Integrator& integrator = default_integrator();
    Dagger dagger(config_.dip_heights, config_.dip_end_times);
    ProductionTracker tracker(config_, *field, integrator, dagger, rng);

    // 6. Set output path
    std::string out_path = config_.output_prefix + "_rank" + std::to_string(rank) + ".csv";
    CsvResultWriter writer(out_path);

    // 7. Simulation loop
    size_t completed = 0;
    const int bar_width = 50;
    const size_t update_interval = std::max(size_t(1), my_ntraj / 20);

    while (source->has_next()) {
        State s = source->next();
        Result res = tracker.run(s);
        writer.write(res);
        completed++;

        if (rank == 0 && (completed % update_interval == 0 || completed == my_ntraj)) {
            float progress = static_cast<float>(completed) / my_ntraj;
            int pos = static_cast<int>(bar_width * progress);

            std::cout << "\r[" ;
            for (int i = 0; i < bar_width; ++i) {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0) << "%" << std::flush;
        }
    }

    if (rank == 0) {
        std::cout << std::endl;
    }

    return 0;
}

void Runner::print_config(SimulationConfig config, int size, int total_traj) const {
    std::cout << "==================================================" << std::endl;
    std::cout << "       UCNtrap Modern Runner Initialized         " << std::endl;
    std::cout << "==================================================" << std::endl;
    // --- 1. Simulation Control ---
    std::cout << " [1. Simulation Control]" << std::endl;
    std::cout << "  - dt:              " << config_.dt << " s" << std::endl;
    std::cout << "  - ntraj:           " << config_.ntraj << std::endl;
    std::cout << "  - seed:            " << config_.seed << std::endl;
    std::cout << "  - MPI Ranks:       " << size << std::endl;

    // --- 2. Simulation Specs ---
    std::cout << "\n [2. Simulation Specs]" << std::endl;
    std::cout << "  - Source:          " << config_.source_type << std::endl;
    std::cout << "  - Field Model:     " << config_.field_model << std::endl;
    std::cout << "  - Integrator:      " << config_.integrator << std::endl;
    std::cout << "  - Tracker Logic:   " << config_.tracker << std::endl;
    
    std::string d_mode = (config_.dagger_mode == DaggerMode::Fast) ? "Fast" :
                            (config_.dagger_mode == DaggerMode::Slow) ? "Slow" : "Segmented";
    std::cout << "  - Dagger Mode:     " << d_mode << std::endl;

    // --- 3. I/O Configuration ---
    std::cout << "\n [3. I/O Configuration]" << std::endl;
    std::cout << "  - Output Prefix:   " << config_.output_prefix << ".csv" << std::endl;
    if (config_.source_type == "pentrack") {
        std::cout << "  - Init File:       " << config_.neutron_init_file << std::endl;
    }

    // --- 4. Experiment Parameters ---
    std::cout << "\n [4. Experiment Parameters]" << std::endl;
    std::cout << "  - Cleaning Time:   " << config_.cleaning_time << " s" << std::endl;
    std::cout << "  - Cleaning H:      " << config_.cleaning_height << " m" << std::endl;
    std::cout << "  - Raised Clean H:  " << config_.raised_cleaning_height << " m" << std::endl;
    std::cout << "  - Hold Time:       " << config_.hold_time << " s" << std::endl;

    // --- 5. Dagger Movement ---
    std::cout << "\n [5. Dagger Movement]" << std::endl;
    std::cout << "  - Dip Heights:     [ ";
    for (double h : config_.dip_heights) std::cout << h << ", ";
    std::cout << "] m" << std::endl;
    std::cout << "  - Dip End Times:   [ ";
    for (double t : config_.dip_end_times) std::cout << t << ", ";
    std::cout << "] s (From hold time)" << std::endl;

    // --- 6. Core Physics Parameters ---
    std::cout << "\n [6. Core Physics Parameters]" << std::endl;
    std::cout << "  - E-Cut:           " << config_.ecut << std::endl;
    std::cout << "  - E-Pow:           " << config_.epow << std::endl;
    std::cout << "  - E-Clean:         " << config_.eclean << " J" << std::endl;
    std::cout << "  - Zeta-Cut:        " << config_.zetacut << " m" << std::endl;
    std::cout << "  - B-Thick:         " << config_.bthick << " nm" << std::endl;
    std::cout << "  - Defect Prob:     " << config_.defect << std::endl;
    std::cout << "  - Heat Mult:       " << config_.heat_mult << std::endl;
    std::cout << "  - Wall Loss Prob:  " << config_.wall_loss_prob << std::endl;
    
    std::cout << "==================================================\n" << std::endl;
    std::cout << "Starting Simulation..." << std::endl;
}

} // namespace ucntrap