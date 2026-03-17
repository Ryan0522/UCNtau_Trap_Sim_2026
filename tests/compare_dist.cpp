#include "ucntrap/runner.hpp"
#include "ucntrap/source/pentrack_reader.hpp"
#include "ucntrap/config.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>

int main() {
    // 1. Setup Config
    ucntrap::SimulationConfig config;
    config.dt = 0.001;
    config.ntraj = 1000;
    
    // 2. Initialize Runner and Source
    auto reader = std::make_unique<ucntrap::PenTrackReader>("neutrons_init.out");
    ucntrap::Runner runner(config, std::move(reader));

    std::cout << "num trace bins: 71657" << std::endl;

    for (int i = 0; i < config.ntraj; ++i) {
        // Run individual neutron
        auto result = runner.run(); 

        // Match your legacy output format:
        // index, t_final, E_final, code, last_x, last_y, last_z
        std::cout << i << ", " 
                  << std::scientific << std::setprecision(6)
                  << result.t << ", " 
                  << result.e_final << ", " 
                  << result.code << ", "
                  << result.x << ", " 
                  << result.y << ", " 
                  << result.z << std::endl;
    }

    return 0;
}