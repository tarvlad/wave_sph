#pragma once

#include <cstddef>
#include "kernel_type.hpp"
#include "misc.hpp"

constexpr double PI = 3.1415926535897931l;


// Request params
constexpr double rho0 = 1.l;
constexpr double delta = 1e-3l;
constexpr double cs = 1.l;
constexpr double extra_border = 75e-2;
constexpr double l_bound = -extra_border;
constexpr double r_bound = 1.l + extra_border;
constexpr size_t n_gas_particles_per_unit = 83;
constexpr double h = 6e-1;
constexpr double time_moment = 2e0;
constexpr double tau = 1e-5L;
constexpr double k = 2e0l * PI;
constexpr double m_g = rho0 / n_gas_particles_per_unit;
constexpr KERNEL_TYPE kernel_type = C6O4;
constexpr size_t n_particles = static_cast<size_t>(
    __c_ceil(n_gas_particles_per_unit * (r_bound - l_bound))
);
constexpr RUNTIME_DENSITY_COMP_SCHEME density_scheme = SCHEME_16_SARANSK;
constexpr RUNTIME_VELOCITY_COMP_SCHEME velocity_scheme = DIRECT_PRESSURE_GRAD_APPROX;
constexpr double phi = 1.l / (h * n_gas_particles_per_unit);

//K = 1 / h
//phi = 1 / (h * n)