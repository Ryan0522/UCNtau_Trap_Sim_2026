#include "ucntrap/config.hpp"
#include <yaml-cpp/yaml.h>
#include <stdexcept>
#include <vector>

namespace ucntrap {

SimulationConfig load_config(const std::string& path) {
    YAML::Node node = YAML::LoadFile(path);
    SimulationConfig config;

    // --- 1. simulation control ---
    config.dt = node["simulation"]["dt"].as<double>();
    config.ntraj = node["simulation"]["ntraj"].as<std::size_t>();
    config.seed = node["simulation"]["seed"].as<unsigned long long>();

    // --- 2. simulation configurations ---
    config.source_type = node["components"]["source_type"].as<std::string>();
    config.field_model = node["components"]["field_model"].as<std::string>();
    config.integrator = node["components"]["integrator"].as<std::string>();
    config.tracker = node["components"]["tracker"].as<std::string>();
    config.writer = node["components"]["writer"].as<std::string>();

    std::string d_mode = node["components"]["dagger_mode"].as<std::string>();
    if (d_mode == "Fast") config.dagger_mode = DaggerMode::Fast;
    else if (d_mode == "Slow") config.dagger_mode = DaggerMode::Slow;
    else if (d_mode == "Segmented") config.dagger_mode = DaggerMode::Segmented;

    // --- 3. file path and IO ---
    config.neutron_init_file = node["io"]["neutron_init_file"].as<std::string>();
    config.x_trace_file = node["io"]["x_trace_file"].as<std::string>();
    config.y_trace_file = node["io"]["y_trace_file"].as<std::string>();
    config.z_trace_file = node["io"]["z_trace_file"].as<std::string>();
    config.output_prefix = node["io"]["output_prefix"].as<std::string>();

    // --- 4. experiment parameters ---
    config.cleaning_time = node["experiment"]["cleaning_time"].as<double>();
    config.cleaning_height = node["experiment"]["cleaning_height"].as<double>();
    config.raised_cleaning_height = node["experiment"]["raised_cleaning_height"].as<double>();
    config.hold_time = node["experiment"]["hold_time"].as<double>();

    // --- 5. dagger movements (relative_dip_end_times is calulated from after hold_time) --- 
    config.dip_heights = node["dagger_motion"]["dip_heights"].as<std::vector<double>>();
    std::vector<double> relative_times = node["dagger_motion"]["rel_dip_end_times"].as<std::vector<double>>();
    config.dip_end_times.clear();
    for (double t_rel : relative_times) {
        config.dip_end_times.push_back(config.hold_time + t_rel);
    }

    // --- 6. physics parameters ---
    config.ecut = node["physics"]["ecut"].as<double>();
    config.epow = node["physics"]["epow"].as<double>();
    config.eclean = node["physics"]["eclean"].as<double>();
    config.thetapow = node["physics"]["thetapow"].as<double>();
    config.zetacut = node["physics"]["zetacut"].as<double>();
    config.bthick = node["physics"]["bthick"].as<double>();
    config.defect = node["physics"]["defect"].as<double>();
    config.heat_mult = node["physics"]["heat_mult"].as<double>();
    config.wall_loss_prob = node["physics"]["wall_loss_prob"].as<double>();

    return config;
}

} // namespace ucntrap