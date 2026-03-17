#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "ucntrap/state.hpp"
#include "ucntrap/numerics/integrator.hpp"
#include "ucntrap/physics/planar_halbach_field.hpp"
#include "ucntrap/physics/trap_halbach_field.hpp"
#include "ucntrap/constants.hpp"
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
        
        if (std::abs(temp[pivot][i]) < 1e-60) return 0.0;

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
    ucntrap::Trace trace = ucntrap::load_trace("./data/xvals.bin", "./data/yvals.bin", "./data/zvals.bin");
    ucntrap::TrapHalbachField field(1.0, trace.x, trace.y, trace.z);
    // ucntrap::PlanarHalbachField field(1.35, 0.05114, 0.0254, 3);
    const auto& integrator = ucntrap::default_integrator();
    
    double dt = 0.00001;
    double eps_pos = 1e-10;
    double eps_mom = 1e-37;

    ucntrap::State base = {0.044, -0.322, -1.138, 1e-27, 2e-28, 5e-29};
    
    std::vector<ucntrap::State> plus_neighbors(6, base);
    std::vector<ucntrap::State> minus_neighbors(6, base);

    plus_neighbors[0].x += eps_pos; plus_neighbors[1].y += eps_pos; plus_neighbors[2].z += eps_pos;
    plus_neighbors[3].px += eps_mom; plus_neighbors[4].py += eps_mom; plus_neighbors[5].pz += eps_mom;

    minus_neighbors[0].x -= eps_pos; minus_neighbors[1].y -= eps_pos; minus_neighbors[2].z -= eps_pos;
    minus_neighbors[3].px -= eps_mom; minus_neighbors[4].py -= eps_mom; minus_neighbors[5].pz -= eps_mom;

    std::cout << "t,det_jacobian" << std::endl;

    double t = 0.0;
    for (int step = 0; step <= 100000; ++step) {
        double J[6][6];
        
        const double q_s = 1.0; 
        const double p_s = ucntrap::constants::kPScale; // 1.25e-27
        double scales[6] = {q_s, q_s, q_s, p_s, p_s, p_s};

        for (int j = 0; j < 6; ++j) {
            double divisor = 2.0 * ((j < 3) ? eps_pos : eps_mom);
            
            double derivatives[6];
            derivatives[0] = (plus_neighbors[j].x - minus_neighbors[j].x) / divisor;
            derivatives[1] = (plus_neighbors[j].y - minus_neighbors[j].y) / divisor;
            derivatives[2] = (plus_neighbors[j].z - minus_neighbors[j].z) / divisor;
            derivatives[3] = (plus_neighbors[j].px - minus_neighbors[j].px) / divisor;
            derivatives[4] = (plus_neighbors[j].py - minus_neighbors[j].py) / divisor;
            derivatives[5] = (plus_neighbors[j].pz - minus_neighbors[j].pz) / divisor;
            for (int i = 0; i < 6; ++i) {
                J[i][j] = derivatives[i] * (scales[j] / scales[i]);
            }
        }

        double det = calculate_determinant(J);

        if (step % 500 == 0) {
            std::cout << std::fixed << std::setprecision(4) << t << "," 
                      << std::setprecision(10) << det << std::endl;
        }

        if (std::abs(det - 1.0) > 1e-3) { 
            std::cerr << "CRITICAL: Liouville violation at t=" << t << " (det=" << det << ")" << std::endl;
            if (t > 0.8) return 1; 
            return 1;
        }

        integrator.step(base, t, dt, field);
        for (int i = 0; i < 6; ++i) {
            integrator.step(plus_neighbors[i], t, dt, field);
            integrator.step(minus_neighbors[i], t, dt, field);
        }
        t += dt;        t += dt;
    }
    
    std::cout << "SUCCESS: Liouville Test Passed." << std::endl;
    return 0;
}