#include "ucntrap/io/trace_loader.hpp"

#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ucntrap {
    
namespace {

std::vector<double> read_raw_double_file(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input.is_open()) {
        throw std::runtime_error("load_trace: unable to open file '" + path + "'");
    }

    std::vector<double> values;
    double value = 0.0;
    while (input.read(reinterpret_cast<char*>(&value), sizeof(double))) {
        values.push_back(value);
    }

    if (!input.eof()) {
        throw std::runtime_error("load_trace: failed while reading file '" + path + "'");
    }

    return values;
}

} // namespace

Trace load_trace(const std::string& x_file,
                 const std::string& y_file,
                 const std::string& z_file) {
    Trace trace;
    trace.x = read_raw_double_file(x_file);
    trace.y = read_raw_double_file(y_file);
    trace.z = read_raw_double_file(z_file);

    if (trace.x.size() != trace.y.size() || trace.x.size() != trace.z.size()) {
        throw std::runtime_error("load_trace: sample length mismatch across X/Y/Z trace files");
    }

    if (trace.x.empty()) {
        throw std::runtime_error("load_trace: trace is empty");
    }

    return trace;
}

} // namespace ucntrap 