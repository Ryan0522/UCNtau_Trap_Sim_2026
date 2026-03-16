#include "ucntrap/source/pentrack_reader.hpp"
#include "ucntrap/constants.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ucntrap {
    
namespace {

std::string resolve_filename(const SimulationConfig& config) {
    if (!config.neutron_init_file.empty()) {
        return config.neutron_init_file;
    }
    return "neutrons_init.out";
}

} // namespace

PenTrackReader::PenTrackReader(const SimulationConfig& config)
    : PenTrackReader(resolve_filename(config), config.ntraj, 1) {}

PenTrackReader::PenTrackReader(std::string filename,
                                std::size_t ntraj,
                                std::size_t skip_lines) {
    std::ifstream input(filename);
    if (!input.is_open()) {
        throw std::runtime_error(
            "PenTrackReader: unable to open file '" + filename + "'"
        );
    }

    std::string line;
    for (std::size_t i = 0; i < skip_lines; ++i) {
        if (!std::getline(input, line)) {
            throw std::runtime_error(
                "PenTrackReader: file ended before header skip completed"
            );
        }
    }

    states_.reserve(ntraj);
    for (std::size_t i = 0; i < ntraj; ++i) {
        if (!std::getline(input, line)) {
            throw std::runtime_error(
                "PenTrackReader: not enough neutrons in input file");
        }
        states_.push_back(parse_line(line));
    }
}

bool PenTrackReader::has_next() const {
    return index_ < states_.size();
}

State PenTrackReader::next() {
    if (!has_next()) {
        throw std::out_of_range("PenTrackReader::next() called after exhaustion");
    }
    return states_[index_++];
}

State PenTrackReader::parse_line(const std::string& line) {
    std::istringstream iss(line);
    std::vector<double> values;
    values.reserve(6);

    double value = 0.0;
    while (iss >> value) {
        values.push_back(value);
    }

    if (values.size() != 6) {
        throw std::runtime_error(
            "PenTrackReader: expected 6 whitespace-separated values per line");
    }

    State s{};
    s.x  = -values[0];
    s.y  =  values[1];
    s.z  =  values[2];
    s.px = -values[3] * constants::kMassN;
    s.py =  values[4] * constants::kMassN;
    s.pz =  values[5] * constants::kMassN;
    return s;
}

} // namespace ucntrap 