#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <mpi.h>
#include <vector>
#include <string>
#include "ucntrap/config.hpp"
#include "ucntrap/runner.hpp"

namespace py = pybind11;
using namespace ucntrap;

void init_mpi_with_args(std::vector<std::string> args) {
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized) {
        int argc = static_cast<int>(args.size());
        
        std::vector<char*> argv_data;
        for (auto& s : args) {
            argv_data.push_back(const_cast<char*>(s.c_str()));
        }
        argv_data.push_back(nullptr);

        char** argv = argv_data.data();
        MPI_Init(&argc, &argv);
    }
}

void finalize_mpi() {
    int initialized, finalized;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    if (initialized && !finalized) {
        MPI_Finalize();
    }
}

PYBIND11_MODULE(ucntrap_py, m) {
    m.def("init_mpi", &init_mpi_with_args, "Initialize MPI with arguments",
          py::arg("args") = std::vector<std::string>());
    m.def("finalize_mpi", &finalize_mpi, "Finalize MPI");

    m.doc() = "UCNtau Simulation Python Bindings";

    // Converting C++ cout to Python sys.stdout
    py::add_ostream_redirect(m, "ostream_redirect");

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
        .def_readwrite("heat_mult", &SimulationConfig::heat_mult)
        .def_readwrite("wall_loss_prob", &SimulationConfig::wall_loss_prob)
        .def_readwrite("cleaning_time", &SimulationConfig::cleaning_time)
        .def_readwrite("cleaning_height", &SimulationConfig::cleaning_height)
        .def_readwrite("raised_cleaning_height", &SimulationConfig::raised_cleaning_height)
        .def_readwrite("holding_time", &SimulationConfig::hold_time)
        .def_readwrite("dip_end_times", &SimulationConfig::dip_end_times)
        .def_readwrite("dagger_mode", &SimulationConfig::dagger_mode)
        .def_readwrite("neutron_init_file", &SimulationConfig::neutron_init_file)
        .def_readwrite("x_trace_file", &SimulationConfig::x_trace_file)
        .def_readwrite("y_trace_file", &SimulationConfig::y_trace_file)
        .def_readwrite("z_trace_file", &SimulationConfig::z_trace_file)
        .def_readwrite("output_prefix", &SimulationConfig::output_prefix)
        .def_readwrite("array_offset", &SimulationConfig::array_offset);

    // Runner
    py::class_<Runner>(m, "Runner")
        .def(py::init<SimulationConfig>())
        .def("run", &Runner::run);
}