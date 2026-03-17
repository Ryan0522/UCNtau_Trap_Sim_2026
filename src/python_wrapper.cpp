#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "ucntrap/config.hpp"
#include "ucntrap/runner.hpp"

namespace py = pybind11;
using namespace ucntrap;

PYBIND11_MODULE(ucntrap_py, m) {
    m.doc() = "UCNtau Simulation Python Bindings";

    // DaggerMode
    py::enum_<DaggerMode>(m, "DaggerMode")
        .value("Fast", DaggerMode::Fast)
        .value("Slow", DaggerMode::Slow)
        .value("Segmented", DaggerMode::Segmented);

    // SimulationConfig
    py::class_<SimulationConfig>(m, "SimulationConfig")
        .def(py::init<>())
        .def_readwrite("dt", &SimulationConfig::dt)
        .def_readwrite("ntraj", &SimulationConfig::ntraj)
        .def_readwrite("seed", &SimulationConfig::seed)
        .def_readwrite("defect", &SimulationConfig::defect)
        .def_readwrite("wall_loss_prob", &SimulationConfig::wall_loss_prob)
        .def_readwrite("cleaning_time", &SimulationConfig::cleaning_time)
        .def_readwrite("holding_time", &SimulationConfig::hold_time)
        .def_readwrite("dagger_mode", &SimulationConfig::dagger_mode)
        .def_readwrite("output_prefix", &SimulationConfig::output_prefix);

    // Runner
    py::class_<Runner>(m, "Runner")
        .def(py::init<SimulationConfig>())
        .def("run", &Runner::run);
}