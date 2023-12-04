#include <iostream>
#include "constants.hpp"
#include "init.hpp"
#include "runtime.hpp"
#include <iomanip>


static_assert(sizeof(int) == 4, "");
int main() {
    std::cout << std::setprecision(16);
    std::cout << std::scientific;
    PARTICLE *particles = new PARTICLE[n_particles];
    init_sph(particles);

    double init_rho_err = 0.l;
    for (size_t i = 0; i < n_particles; i++) {
        PARTICLE &p = particles[i];
        if (0.l <= p.coordinate && p.coordinate <= 1.l) {
            double error = std::abs(p.density - rho_analytical(p.coordinate, 0.l));
            if (error > init_rho_err) {
                init_rho_err = error;
            }
        }
    }

    constexpr size_t n_time = __c_ceil(time_moment / tau);
    for (size_t t = 0; t < n_time; t++) {
        right_side_update(particles);
        update(particles);

        if (OUT_CONDITION(t)) {
            CLR_SCREEN;
            std::cout << "Step " << t << " of " << n_time << std::endl;
        }
    }

    double max_error_v = 0.l;
    double max_error_rho = 0.l;
    for (size_t i = 0; i < n_particles; i++) {
        PARTICLE &p = particles[i];
        if (p.coordinate >= 0.l && p.coordinate <= 1.l) {
            double v_err = std::abs(p.velocity - v_analytical(p.coordinate, time_moment));
            if (v_err > max_error_v) {
                max_error_v = v_err;
            }

            double rho_err = std::abs(p.density - rho_analytical(p.coordinate, time_moment));
            if (rho_err > max_error_rho) {
                max_error_rho = rho_err;
            }
        }
    }

    CLR_SCREEN;
    std::cout << "phi = " << phi << std::endl;
    std::cout << "max_error_v = " << max_error_v << std::endl;
    std::cout << "max_error_rho = " << max_error_rho << std::endl;

    delete[] particles;
    return 0;
}
