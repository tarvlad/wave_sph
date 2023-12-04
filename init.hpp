#pragma once

#include <iostream>
#include "constants.hpp"
#include "analytic.hpp"
#include "particle.hpp"


inline
void density_based_arrangement(PARTICLE *particles) {
    double *borders = new double[n_particles + 1];
    borders[0] = l_bound;
    borders[n_particles] = r_bound;

    //TODO tau or pow(10, 5) ??? 
    double delta_b = (r_bound - l_bound) * tau / n_particles;

    CLR_SCREEN;
    std::cout << "Border " << 0 << " of " << n_particles << " done" << std::endl;
    for (size_t i = 1; i < n_particles; i++) {
        double curr_b = borders[i - 1];
        double curr_i = 0.l;
        do {
            curr_b += delta_b;
            double coord = curr_b;
            curr_i += delta_b * rho_analytical(coord, 0.l);
        } while (curr_i < m_g);
        borders[i] = curr_b;

        if (OUT_CONDITION(i)) {
            CLR_SCREEN;
            std::cout << "Border " << i << " of " << n_particles << " done" << std::endl;
        }
    }
    CLR_SCREEN;
    std::cout << "Border " << n_particles << " of " << n_particles << " done" << std::endl;

    for (size_t i = 0; i < n_particles; i++) {
        double coord = (borders[i] + borders[i + 1]) / 2.l;
        particles[i].coordinate = coord;
    }

    delete[] borders;
}


inline
void init_sph(PARTICLE *particles) {
    density_based_arrangement(particles);

    for (size_t i = 0; i < n_particles; i++) {
        PARTICLE &p = particles[i];
        p.velocity = v_analytical(p.coordinate, 0.l);
        p.density = rho_analytical(p.coordinate, 0.l);

        p.r_coordinate = 0.l;
        p.r_velocity = 0.l;
        p.r_density = 0l;
    }
}