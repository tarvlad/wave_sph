#pragma once

#include "constants.hpp"
#include <cmath>

inline
double v_analytical(double x, double t) {
    return delta * (cs / rho0) * cos(k * (x - cs * t));
}


inline
double rho_analytical(double x, double t) {
    return rho0 + delta * cos(k * (x - cs * t));
}