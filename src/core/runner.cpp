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
        std::cout << "Starting Simulation with " << size << " ranks." << std::endl;
        std::cout << "Total Trajectories: " << total_ntraj << std::endl;
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
    ProductionTracker tracker(config_, *field, integrator, dagger);

    // 6. Set output path
    std::string out_path = config_.output_prefix + "_rank" + std::to_string(rank) + ".csv";
    CsvResultWriter writer(out_path);

    // 7. Simulation loop
    while (source->has_next()) {
        State s = source->next();
        Result res = tracker.run(s);
        writer.write(res);
    }

    return 0;
}

} // namespace ucntrap