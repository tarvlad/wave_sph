#include "Vector_extension.hpp"
#include <iostream>
#include <cmath>

std::vector<double> operator+(const std::vector<double> &a, const std::vector<double> &b) {
    if (a.size() != b.size()) {
        std::cout << "a.size() != b.size()" << std::endl;
        exit(1);
    }
    std::vector<double> c(a.size());
    for (int i = 0; i < a.size(); ++i) {
        c[i] = a[i] + b[i];
    }
    return c;
}

std::vector<double> operator-(const std::vector<double> &a, const std::vector<double> &b) {
    if (a.size() != b.size()) {
        std::cout << "a.size() != b.size()" << std::endl;
        exit(1);
    }
    std::vector<double> c(a.size());
    for (int i = 0; i < a.size(); ++i) {
        c[i] = a[i] - b[i];
    }
    return c;
}

double operator*(const std::vector<double> &a, const std::vector<double> &b) {
    if (a.size() != b.size()) {
        std::cout << "a.size() != b.size()" << std::endl;
        exit(1);
    }
    double c = 0;
    for (int i = 0; i < a.size(); ++i) {
        c += a[i] * b[i];
    }
    return c;
}

std::vector<double> operator*(const double &b, const std::vector<double> &a) {
    std::vector<double> c(a.size());
    for (int i = 0; i < a.size(); ++i) {
        c[i] = a[i] * b;
    }
    return c;
}

std::vector<double> operator*(const std::vector<double> &a, const double &b) {
    std::vector<double> c(a.size());
    for (int i = 0; i < a.size(); ++i) {
        c[i] = a[i] * b;
    }
    return c;
}

std::vector<double> operator/(const std::vector<double> &a, const double &b) {
    std::vector<double> c(a.size());
    for (int i = 0; i < a.size(); ++i) {
        c[i] = a[i] / b;
    }
    return c;
}

double abs(const std::vector<double> &a) {
    double y = 0;
    y = sqrt(a * a);
    return y;
}
