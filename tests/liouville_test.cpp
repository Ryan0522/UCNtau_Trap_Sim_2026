#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "ucntrap/state.hpp"
#include "ucntrap/numerics/integrator.hpp"
#include "ucntrap/physics/planar_halbach_field.hpp"

// Manual 6x6 Determinant using Gaussian Elimination
double calculate_determinant(double mat[6][6]) {
    double det = 1.0;
    double temp[6][6];
    for(int i=0; i<6; ++i) for(int j=0; j<6; ++j) temp[i][j] = mat[i][j];

    for (int i = 0; i < 6; ++i) {
        int pivot = i;
        for (int j = i + 1; j < 6; ++j) {
            if (std::abs(temp[j][i]) > std::abs(temp[pivot][i])) pivot = j;
        }
        if (std::abs(temp[pivot][i]) < 1e-20) return 0.0;

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
    // 1. Setup Environment (Using Planar field for a clean benchmark)
    // Parameters: B_rem, space, thick, terms
    ucntrap::PlanarHalbachField field(1.35, 0.05114, 0.0254, 3);
    const auto& integrator = ucntrap::default_integrator();
    
    double dt = 0.001;
    double epsilon = 1e-8; // Perturbation size for Jacobian

    // 2. Initial Base State (x, y, z, px, py, pz)
    ucntrap::State base = {0.044, -0.322, -1.138, 1e-27, 2e-28, 5e-29};
    
    // Create 6 perturbed neighbors representing the 6 phase-space axes
    std::vector<ucntrap::State> neighbors(6, base);
    neighbors[0].x += epsilon;
    neighbors[1].y += epsilon;
    neighbors[2].z += epsilon;
    neighbors[3].px += epsilon;
    neighbors[4].py += epsilon;
    neighbors[5].pz += epsilon;

    std::cout << "t,det_jacobian" << std::endl;

    double t = 0.0;
    for (int step = 0; step <= 5000; ++step) {
        // Construct the Jacobian Matrix: J_ij = (pos_i_neighbor_j - pos_i_base) / epsilon
        double J[6][6];
        for (int j = 0; j < 6; ++j) {
            J[0][j] = (neighbors[j].x - base.x) / epsilon;
            J[1][j] = (neighbors[j].y - base.y) / epsilon;
            J[2][j] = (neighbors[j].z - base.z) / epsilon;
            J[3][j] = (neighbors[j].px - base.px) / epsilon;
            J[4][j] = (neighbors[j].py - base.py) / epsilon;
            J[5][j] = (neighbors[j].pz - base.pz) / epsilon;
        }

        double det = calculate_determinant(J);
        
        if (step % 500 == 0) {
            std::cout << std::fixed << std::setprecision(3) << t << "," 
                      << std::setprecision(10) << det << std::endl;
        }

        // Check for major drift (Symplectic integrators should stay near 1.0)
        if (std::abs(det - 1.0) > 1e-4) {
            std::cerr << "CRITICAL: Liouville violation at t=" << t << " (det=" << det << ")" << std::endl;
            return 1;
        }

        // Evolve the system
        integrator.step(base, t, dt, field);
        for (auto& n : neighbors) integrator.step(n, t, dt, field);
        t += dt;
    }

    std::cout << "SUCCESS: Phase-space volume preserved over 5 seconds." << std::endl;
    return 0;
}