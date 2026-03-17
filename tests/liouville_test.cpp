#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "ucntrap/state.hpp"
#include "ucntrap/numerics/integrator.hpp"
#include "ucntrap/physics/planar_halbach_field.hpp"
#include "ucntrap/physics/trap_halbach_field.hpp"
#include "ucntrap/io/trace_loader.hpp"

double calculate_determinant(double mat[6][6]) {
    double det = 1.0;
    double temp[6][6];
    for(int i=0; i<6; ++i) for(int j=0; j<6; ++j) temp[i][j] = mat[i][j];

    for (int i = 0; i < 6; ++i) {
        int pivot = i;
        for (int j = i + 1; j < 6; ++j) {
            if (std::abs(temp[j][i]) > std::abs(temp[pivot][i])) pivot = j;
        }
        if (std::abs(temp[pivot][i]) < 1e-22) return 0.0;

        if (pivot != i) {
            for (int k = 0; k < 6; ++k) std::swap(temp[i][k], temp[pivot][k]);
            det *= -1.0;
        }

        det *= temp[i][i];
        for (int j = i + 1; j < 6; ++j) {
            double factor = temp[j][i] / temp[i][i];
            for (int k = i + 1; k < 6; ++k) temp[j][k] -= factor * temp[i][k];
        }
    }
    return det;
}

int main() {
    ucntrap::PlanarHalbachField field(1.35, 0.05114, 0.0254, 3);
    // auto trace = ucntrap::load_trace("./data/xvals.bin", "./data/yvals.bin", "./data/zvals.bin");
    // ucntrap::TrapHalbachField field(1.0, trace.x, trace.y, trace.z);
    const auto& integrator = ucntrap::default_integrator();
    
    double dt = 0.0001;
    double eps_pos = 1e-8; // perturbation
    double eps_mom = 1e-35;

    ucntrap::State base = {0.044, -0.322, -1.138, 1e-27, 2e-28, 5e-29};
    
    std::vector<ucntrap::State> neighbors(6, base);
    neighbors[0].x += eps_pos; neighbors[1].y += eps_pos; neighbors[2].z += eps_pos;
    neighbors[3].px += eps_mom; neighbors[4].py += eps_mom; neighbors[5].pz += eps_mom;

    std::cout << "t,det_jacobian" << std::endl;

    double t = 0.0;
    for (int step = 0; step <= 50000; ++step) {
        double J[6][6];
        for (int j = 0; j < 6; ++j) {
            double divisor = (j < 3) ? eps_pos : eps_mom;
            J[0][j] = (neighbors[j].x - base.x) / divisor;
            J[1][j] = (neighbors[j].y - base.y) / divisor;
            J[2][j] = (neighbors[j].z - base.z) / divisor;
            J[3][j] = (neighbors[j].px - base.px) / divisor;
            J[4][j] = (neighbors[j].py - base.py) / divisor;
            J[5][j] = (neighbors[j].pz - base.pz) / divisor;
        }

        double det = calculate_determinant(J);
        
        if (step % 500 == 0) {
            std::cout << std::fixed << std::setprecision(3) << t << "," 
                      << std::setprecision(10) << det << std::endl;
        }

        if (std::abs(det - 1.0) > 1e-1) {
            std::cerr << "CRITICAL: Liouville violation at t=" << t << " (det=" << det << ")" << std::endl;
            return 1;
        }

        integrator.step(base, t, dt, field);
        for (auto& n : neighbors) integrator.step(n, t, dt, field);
        t += dt;
    }
    std::cout << "SUCCESS: Liouville Test Passed." << std::endl;
    return 0;
}