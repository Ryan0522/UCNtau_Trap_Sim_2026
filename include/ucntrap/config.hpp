#pragma once
#include <cstddef>
#include <string>
#include <vector>

namespace ucntrap {
 
struct SimulationConfig {
    double dt = 1.0e-4;
    std::size_t ntraj = 1000;
    unsigned long long seed = 1;

    std::string source_type = "pentrack";
    std::string field_model = "trap";
    std::string integrator = "symplectic";
    std::string tracker = "production"; // production, cleaning, Lyap (legacy) 
    std::string writer = "legacy_binary";

    std::string neutron_init_file;
    std::string x_trace_file;
    std::string y_trace_file;
    std::string z_trace_file;
    std::string output_prefix = "results/run";

    double cleaning_time = 50.0;
    double cleaning_height = 0.38;
    double raised_cleaning_height = 0.43;
    double hold_time = 20.0; // Calculate First dip time from here

    std::vector<double> dip_heights;
    std::vector<double> dip_end_times;

    double ecut = 7.2092;
    double epow = 1.17727;
    double eclean = 5.571749397933261e-27;
    double thetapow = 0.275457;
    double zetacut = 0.0127203;
    double bthick = 5.59909;
    double defect = 7.0e-5;
    double heat_mult = 1.0;
};

} // namespace ucntrap