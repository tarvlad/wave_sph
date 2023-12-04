#pragma once


#include "constants.hpp"
#include <cmath>


inline
double w(double x_a, double x_b) {
    double q = std::abs(x_a - x_b) / h;

    double y = 0;
    if (q > 1.l) {
        return y;
    }

    if constexpr (kernel_type == C4O2) {
        y = 3.l / (2.l * h);
        y *= pow(1.l - q, 5);
        y *= 8.l * pow(q, 2) + 5.l * q + 1.l;
    }
    if constexpr (kernel_type == C4O4) {
        y = 3.l / (4.l * h);
        y *= (60.l - 330.l * pow(q, 2)) / 19.l;
        y *= pow(1.l - q, 5);
        y *= 8.l * pow(q, 2) + 5.l * q + 1;
    }
    if constexpr (kernel_type == C6O2) {
        y = 55.l / (32.l * h);
        y *= pow(1.l - q, 7);
        y *= 1.l + 7.l * q + 19.l * pow(q, 2) + 21.l * pow(q, 3);
    }
    if constexpr (kernel_type == C6O4) {
        y = 429.l / (160.l * h);
        y *= pow(q - 1.l, 7);
        y *= 7.l * pow(q, 2) - 1.l;
        y *= 21.l * pow(q, 3) + 19.l * pow(q, 2) + 7.l * q + 1.l;
    }

    return y;
}


inline
double dw_dq(double x_a, double x_b) {
    double q = std::abs(x_a - x_b) / h;

    double y = 0;
    if (q > 1.l) {
        return y;
    }

    if constexpr (kernel_type == C4O2) {
        y = 3.l / (2.l * h);
        y *= pow(1.l - q, 4);
        y *= -56.l * pow(q, 2) - 14.l * q;
    }
    if constexpr (kernel_type == C4O4) {
        y = 3.l * 30.l * q / (2.l * h * 19.l);
        y *= pow(1.l - q, 4);
        y *= 396.l * pow(q, 3) + 44.l * pow(q, 2) - 100.l * q - 25.l;
    }
    if constexpr (kernel_type == C6O2) {
        y = -165.l / (16.l * h);
        y *= pow(-1.l + q, 6);
        y *= 9;
        y *= (3.l + q * (18.l + 35.l * q));
    }
    if constexpr (kernel_type == C6O4) {
        y = 429.l / (40.l * h);
        y *= pow(q - 1.l, 6);
        y *= q;
        y *= -8.l + q * (-48.l + 7.l * q * (-9.l + q * (26.l + 63.l * q)));
    }

    return y;
}


inline
double grad_w(double x_a, double x_b) {
    double r_ij = x_a - x_b;

    double grad = 0.l;
    if (r_ij != 0.l) {
        grad = r_ij / (std::abs(r_ij) * h) * dw_dq(x_a, x_b);
    }
    return grad;
}