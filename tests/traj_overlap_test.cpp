#include "ucntrap/state.hpp"
#include "ucntrap/numerics/integrator.hpp"
#include "ucntrap/physics/trap_halbach_field.hpp"
#include "ucntrap/physics/planar_halbach_field.hpp"
#include "ucntrap/source/pentrack_reader.hpp"
#include "ucntrap/io/trace_loader.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

void run_overlap_test(const std::string& model_name, const std::string& ref_file_path) {
    // 1. Setup physics environment
    auto trace = ucntrap::load_trace("./data/xvals.bin", "./data/yvals.bin", "./data/zvals.bin");
    const auto& integrator = ucntrap::default_integrator();
    
    std::unique_ptr<ucntrap::FieldModel> field;
    if (model_name == "trap") {
        field = std::make_unique<ucntrap::TrapHalbachField>(1.0, trace.x, trace.y, trace.z);
    } else {
        field = std::make_unique<ucntrap::PlanarHalbachField>(1.35, 0.05114, 0.0254, 3);
    }

    // 2. Initial state
    std::string init_line = "0.0442715218949471 -0.32237448462642004 -1.1380754512554265 -0.733273860557217 0.117996751565768 0.0029265797433262";
    ucntrap::State s = ucntrap::PenTrackReader::parse_line(init_line);

    std::ifstream ref_file(ref_file_path);
    std::ofstream out("./tests/" + model_name + "_comparison.csv");
    out << "t,modern_x,modern_y,modern_z,legacy_x,legacy_y,legacy_z\n";

    std::string header; std::getline(ref_file, header);
    double t = 0.0, dt = 0.001;
    double r_t, r_x, r_y, r_z, r_px, r_py, r_pz;

    while (t <= 50.0 && ref_file >> r_t >> r_x >> r_y >> r_z >> r_px >> r_py >> r_pz) {
        out << t << "," << s.x << "," << s.y << "," << s.z << "," 
            << r_x << "," << r_y << "," << r_z << "\n";
        integrator.step(s, t, dt, *field);
        t += dt;
    }
    std::cout << "Completed: " << model_name << " analysis." << std::endl;
}

int main() {
    run_overlap_test("trap", "./tests/legacy_traj.txt");
    run_overlap_test("planar", "./tests/legacy_traj_planar.txt");
    return 0;
}