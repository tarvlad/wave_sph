#pragma once
#include "constants.hpp"
#include "kernels.hpp"
#include "particle.hpp"

inline
void right_side_update(PARTICLE *particles) {
    for (size_t i = 0; i < n_particles; i++) {
        PARTICLE &p = particles[i];
        p.r_coordinate = p.coordinate + tau * p.velocity;

        if constexpr (density_scheme == SCHEME_16_SARANSK) {
            p.r_density = 0.l;
            for (size_t j = 0; j < n_particles; j++) {
                PARTICLE &p1 = particles[j];
                double r_ij = p.coordinate - p1.coordinate;
                if (std::abs(r_ij) <= h) {
                    p.r_density += (p.velocity + p1.velocity) * grad_w(p.coordinate, p1.coordinate);
                }
            }
            p.r_density *= -m_g;
            p.r_density = p.density + tau * p.r_density;
        }


        if constexpr (velocity_scheme == DIRECT_PRESSURE_GRAD_APPROX) {
            p.r_velocity = 0.l;
            for (size_t j = 0; j < n_particles; j++) {
                PARTICLE &p1 = particles[j];
                double r_ij = p.coordinate - p1.coordinate;
                if (std::abs(r_ij) <= h) {
                    p.r_velocity = p.r_velocity + grad_w(p.coordinate, p1.coordinate);
                }
            }
            p.r_velocity *= (-1.l) * pow(cs, 2) * m_g / p.density;
            p.r_velocity = p.velocity + tau * p.r_velocity;
        }
    }
}

inline
void update(PARTICLE *particles) {
    for (size_t i = 0; i < n_particles; i++) {
        PARTICLE &p = particles[i];
        p.coordinate = p.r_coordinate;
        p.velocity = p.r_velocity;
        p.density = p.r_density;
    }
}