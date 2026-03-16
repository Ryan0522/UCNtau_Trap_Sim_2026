#include "ucntrap/runner.hpp"
#include "ucntrap/config.hpp"

#include <iostream>

namespace ucntrap {

Runner::Runner(SimulationConfig config) : config_(config) {}

int Runner::run() const {
    std::cout << "UCNTrap simulation\n";
    
    std::cout << " ============= Config ============= \n";
    std::cout << "source = " << config_.source_type << "\n";
    std::cout << "field model = " << config_.field_model << "\n";
    std::cout << "integrator = " << config_.integrator << "\n";
    std::cout << "tracker = " << config_.tracker << "\n";
    std::cout << "writer = " << config_.writer << "\n";
    std::cout << "output_prefix = " << config_.output_prefix << "\n";

    std::cout << " ============= Parameters ============= \n";
    std::cout << "dt = " << config_.dt << "\n";
    std::cout << "ntraj = " << config_.ntraj << "\n";
    std::cout << "clean time = " << config_.cleaning_time << "\n";
    std::cout << "clean height = " << config_.cleaning_height << "\n";
    std::cout << "raised clean height = " << config_.raised_cleaning_height << "\n";
    std::cout << "hold time = " << config_.hold_time << "\n";

    std::cout << "dip heights = <";
    for (const auto& i : config_.dip_heights) {
        std::cout << i << ", ";
    }
    std::cout << ">\n";

    std::cout << "dip end times = <";
    for (const auto& i : config_.dip_end_times) {
        std::cout << i << ", ";
    }
    std::cout << ">\n";

    return 0;
}

} // namespace ucntrap