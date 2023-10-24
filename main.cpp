#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip> // для использования setprecision()
#include "Vector_extension.hpp"


#define CLR_SCREEN system("cls")
#define OUTPUT_CONDITION(var) (var % 50 == 0) 

const double PI = acos(-1.0);

//папка и файлы для записи результатов
std::string folder = "C:\\Users\\tarvlad\\source\\WaveSPH";
std::string file_1D = folder + "/results_1D.txt";

// начальные данные
struct initial_data {
//    const double nu = 1.0 / 80;
    const double rho_0 = 1;
    const double Delta = pow(10, -3);
    const double Cs = 1;
//    const double gamma = 7;
    const double L_boundary = -0.7;
    const double R_boundary = 1.7;
    const int N_gas_particles_per_unit_side = 63;
    const double h = 0.8;
    const double moment_in_time = 0.002;
    const double tau = 0.00001;
//    const double eta_squared = 0.001 * pow(h, 2);
    const double k = 2 * PI;
    double m_g = rho_0 / N_gas_particles_per_unit_side;
};

// аналитические функции
double v_x_analytical(const std::vector<double> &X, const double &t, initial_data &InitData);

double rho_analytical(const std::vector<double> &X, const double &t, initial_data &InitData);

//double p_analytical(const std::vector<double> &X, const double &t, initial_data &InitData);

// функции ядра
double W (const std::vector<double> &x, const std::vector<double> &x0, initial_data &InitData);

double dW_dq (const std::vector<double> &x, const std::vector<double> &x0, initial_data &InitData);

//double d2W_dq2 (const std::vector<double> &x, const std::vector<double> &x0, initial_data &InitData);

std::vector<double> grad_W (const std::vector<double> &x, const std::vector<double> &x0, initial_data &InitData);

// частицы
struct particle {
    // искомые величины и правые части для их расчёта
    std::vector<double> coordinate, right_side_coordinate;
    std::vector<double> velocity, right_side_velocity;
    double density, right_side_density;
//    double pressure, right_side_pressure;
};

// функции частиц
void uniform_arrangement(std::vector<particle> &particles, initial_data &InitData);

void density_based_arrangement(std::vector<particle> &particles, initial_data &InitData);

void init_SPH(std::vector<particle> &particles, initial_data &InitData);

void right_side_calculation(std::vector<particle> &particles, initial_data &InitData);

void update(std::vector<particle> &particles);


int main(int argc, const char * argv[]) {
    std::cout << std::boolalpha; //чтобы cout выводил true и false, вместо 1 и 0
    std::cout << std::setprecision(16);

    // вводим начальные данные
    initial_data InitData;
    
    // расчёт вспомогательных констант
    int N_g = InitData.N_gas_particles_per_unit_side * (InitData.R_boundary - InitData.L_boundary);
    
    // вводим частицы газа/
    std::vector<particle> particles(N_g);
    
    // расчитываем начальные значения искомых функций
    init_SPH(particles, InitData);
    
    //считаем максимальную погрешность плотности
    double init_max_error_rho = 0;
    double error = 0;
    for (auto &pi : particles) {
        if (0 <= pi.coordinate[0] && pi.coordinate[0] <= 1) {
            error = abs(pi.density - rho_analytical(pi.coordinate, 0, InitData));
            if (error > init_max_error_rho) {
                init_max_error_rho = error;
            }
        }
    }
    
    // количество шагов по времени
    int N_time = ceil(InitData.moment_in_time / InitData.tau);
    
    // расчитываем изменения функций с течением времени
    for (int n = 1; n <= N_time; n++) {
        right_side_calculation(particles, InitData);
        
        update(particles);
        
        if (OUTPUT_CONDITION(n)) {
            CLR_SCREEN;
            std::cout << "Step " << n << " of " << N_time << std::endl;
        }
    }
    
    //считаем максимальные погрешности скорости и плотности
    double max_error_v = 0, max_error_rho = 0;
    for (auto &pi : particles) {
        if (0 <= pi.coordinate[0] && pi.coordinate[0] <= 1) {
            error = abs(pi.velocity[0] - v_x_analytical(pi.coordinate, InitData.moment_in_time, InitData));
            if (error > max_error_v) {
                max_error_v = error;
            }
            
            error = abs(pi.density - rho_analytical(pi.coordinate, InitData.moment_in_time, InitData));
            if (error > max_error_rho) {
                max_error_rho = error;
            }
        }
    }
    
    // вывод результатов расчёта
    CLR_SCREEN;
    std::cout << "tau = " << InitData.tau << ", h = " << InitData.h << ", N = " << InitData.N_gas_particles_per_unit_side << ", N_time = " << N_time << std::endl;
    std::cout << "End max speed error = " << max_error_v << std::endl;
    std::cout << "Begin max density error = " << init_max_error_rho << std::endl;
    std::cout << "End max density error = " << max_error_rho << std::endl;
    
    //запись результатов в файл
    std::ofstream out_1D;
    out_1D.open(file_1D);
    if (out_1D.is_open()) {
        out_1D << std::setprecision(16);
        out_1D << "x, u, rho, an_u, an_rho" << std::endl;
        for (auto &pi : particles) {
            if (pi.coordinate[0] >= 0.0 && pi.coordinate[0] <= 1.0) {
                out_1D << pi.coordinate[0] << ", " << pi.velocity[0] << ", " << pi.density;
                out_1D << ", " << v_x_analytical(pi.coordinate, InitData.moment_in_time, InitData) << ", " << rho_analytical(pi.coordinate, InitData.moment_in_time, InitData);
                out_1D << std::endl;
            }
        }
        out_1D << std::endl;
    }
    out_1D.close();
    
    std::cout << std::endl << "Finish!" << std::endl;
    return 0;
}

//аналитические функции
double v_x_analytical(const std::vector<double> &X, const double &t, initial_data &InitData) {
    double y = 0;
    y = InitData.Delta * (InitData.Cs / InitData.rho_0) * cos(InitData.k * X[0] - InitData.k * InitData.Cs * t);
    return y;
}

double rho_analytical(const std::vector<double> &X, const double &t, initial_data &InitData) {
    double y = 0;
    y = InitData.rho_0 + InitData.Delta * cos(InitData.k * X[0] - InitData.k * InitData.Cs * t);
    return y;
}

//double p_analytical(const std::vector<double> &X, const double &t, initial_data &InitData) {
//    double y = 0;
//    y = pow(InitData.Cs, 2) * rho_analytical(X, t, InitData);
//    return y;
//}

// функции ядра
double W (const std::vector<double> &x, const std::vector<double> &x0, initial_data &InitData) {
    double h = InitData.h;
    double q = abs(x - x0) / h;
    
    double y = 0;
    if (q <= 1) {
        // S6 O2
//        y = 55 / (32 * h);
//        y *= pow(1 - q, 7);
//        y *= 1 + 7 * q + 19 * pow(q, 2) + 21 * pow(q, 3);
        // S4 O2
        y = 3 / (2 * h);
        y *= pow(1 - q, 5);
        y *= 8 * pow(q, 2) + 5 * q + 1;

        // S4 O4
//        y = 3 / (4 * h);
//        y *= (60 - 330 * pow(q, 2)) / 19;
//        y *= pow(1 - q, 5);
//        y *= 8 * pow(q, 2) + 5 * q + 1;
    }
    return y;
}

double dW_dq (const std::vector<double> &x, const std::vector<double> &x0, initial_data &InitData) {
    double h = InitData.h;
    double q = abs(x - x0) / h;
    
    double y = 0;
    if (q <= 1) {
        // S6 O2
//        y = -165 / (16 * h);
//        y *= pow(-1 + q, 6);
//        y *= q;
//        y *= (3 + q * (18 + 35 * q));
        // S4 O2
        y = 3 / (2 * h);
        y *= pow(1 - q, 4);
        y *= -56 * pow(q, 2) - 14 * q;
        
        // S4 O4
//        y = 3 * 30 * q / (2 * h * 19);
//        y *= pow(1 - q, 4);
//        y *= 396 * pow(q, 3) + 44 * pow(q, 2) - 100 * q - 25;
    }
    return y;
}

//double d2W_dq2 (const std::vector<double> &x, const std::vector<double> &x0, initial_data &InitData) {
//    double h = InitData.h;
//    double q = abs(x - x0) / h;
//
//    double y = 0;
//    if (q <= 1) {
//        // S4 O2
////        y = 3 / (2 * h);
////        y *= pow(1 - q, 3);
////        y *= 336 * pow(q, 2) - 42 * q - 14;
//
//        // S4 O4
//        y = - 3 * 30 / (2 * h * 19);
//        y *= pow(1 - q, 3);
//        y *= 3186 * pow(q, 4) - 1276 * pow(q, 3) - 732 * pow(q, 2) + 75 * q + 25;
//    }
//    return y;
//}

std::vector<double> grad_W (const std::vector<double> &x, const std::vector<double> &x0, initial_data &InitData) {
    double h = InitData.h;
    std::vector<double> r_ij = x - x0;
    
    std::vector<double> c(r_ij.size());
    if (abs(r_ij) != 0.0) {
        c = r_ij / (abs(r_ij) * h) * dW_dq(x, x0, InitData);
    } else {
        c = c * 0.0;
    }
    return c;
}

// функции частиц
void uniform_arrangement(std::vector<particle> &particles, initial_data &InitData) {
    double L = InitData.L_boundary, R = InitData.R_boundary;
    
    unsigned long N = particles.size();
    double step = (R - L) / (N - 1);
    std::vector<double> coord(1);
    coord[0] = L;
    for (int i = 0; i < N; i++) {
        particles[i].coordinate = coord;
        coord[0] += step;
    }
}

void density_based_arrangement(std::vector<particle> &particles, initial_data &InitData) {
    double L = InitData.L_boundary, R = InitData.R_boundary;
    unsigned long N = particles.size();
    std::vector<double> borders(N + 1);
    borders[0] = L;
    borders[N] = R;
    double delta_b = (R - L) / (pow(10, 5) * N);
    std::vector<double> coord(1);
    coord = coord * 0;
    double curr_b;
    double curr_I;
    CLR_SCREEN;
    std::cout << "Border [" << 0 << "] of [" << N << "] done" << std::endl;
    for (int i = 1; i < N; i++) {
        curr_b = borders[i - 1];
        curr_I = 0;
        do {
            curr_b += delta_b;
            coord[0] = curr_b;
            curr_I += delta_b * rho_analytical(coord, 0, InitData);
        } while (curr_I < InitData.m_g);
        borders[i] = curr_b;
        if (OUTPUT_CONDITION(i)) {
            CLR_SCREEN;
            std::cout << "Border [" << i << "] of [" << N << "] done" << std::endl;
        }
    }
    CLR_SCREEN;
    std::cout << "Border [" << N << "] of [" << N << "] done" << std::endl;
    for (int i = 0; i < N; i++) {
        coord[0] = (borders[i] + borders[i + 1]) / 2;
        particles[i].coordinate = coord;
    }
    //std::cout << "Расстановка частиц завершена" << std::endl;
}

void init_SPH(std::vector<particle> &particles, initial_data &InitData) {
    // равномерная расстановка частиц
//    uniform_arrangement(particles, InitData);
    
    // расстановка частиц в соответствии с плотностью
    density_based_arrangement(particles, InitData);
    
    std::vector<double> vel(1);
    for (auto &pi : particles) {
        // аналитическая скорость
        vel[0] = v_x_analytical(pi.coordinate, 0, InitData);
        pi.velocity = vel;
        
        // аналитическая плотность
        pi.density = rho_analytical(pi.coordinate, 0, InitData);
        // численная плотность
//        pi.density = 0;
//        for (auto &pj : particles) {
//            double q = abs(pi.coordinate - pj.coordinate) / InitData.h;
//            if(q <= 1) {
//                pi.density += W(pi.coordinate, pj.coordinate, InitData);
//            }
//        }
//        pi.density *= InitData.m_g;
        
        // аналитическое давление
//        pi.pressure = p_analytical(pi.coordinate, 0, InitData);
        
        // правые части
        pi.right_side_coordinate = pi.coordinate * 0;
        pi.right_side_velocity = pi.velocity * 0;
        pi.right_side_density = pi.density * 0;
//        pi.right_side_pressure = pi.pressure * 0;
    }
}

void right_side_calculation(std::vector<particle> &particles, initial_data &InitData) {
    // расстояния между частицами
    std::vector<double> r_ij(1);
    
    // расчёт правых частей
    for (auto &pi : particles) {
        // координаты
        pi.right_side_coordinate = pi.velocity;
        pi.right_side_coordinate = pi.coordinate + InitData.tau * pi.right_side_coordinate;
        
        // плотности (по фактическому расположению частиц)
//        pi.right_side_density = 0;
//        for (auto &pj : particles) {
//            double q = abs(pi.coordinate - pj.coordinate) / InitData.h;
//            if (q <= 1) {
//                pi.right_side_density += W(pi.coordinate, pj.coordinate, InitData);
//            }
//        }
//        pi.right_side_density *= InitData.m_g;
        
        // плотности (по уравнению неразрывности)
//        pi.right_side_density = 0;
//        for (auto &pj : particles) {
//            r_ij = pi.coordinate - pj.coordinate;
//            if (abs(r_ij) <= InitData.h) {
//                pi.right_side_density += pj.velocity / pj.density * grad_W(pi.coordinate, pj.coordinate, InitData);
//            }
//        }
//        pi.right_side_density *= -pi.density * InitData.m_g;
//        pi.right_side_density = pi.density + InitData.tau * pi.right_side_density;
        
        // плотности (формула 16 из конференции Саранск)
        pi.right_side_density = 0;
        for (auto &pj : particles) {
            r_ij = pi.coordinate - pj.coordinate;
            if (abs(r_ij) <= InitData.h) {
                pi.right_side_density += (pi.velocity + pj.velocity) * grad_W(pi.coordinate, pj.coordinate, InitData);
            }
        }
        pi.right_side_density *= -InitData.m_g;
        pi.right_side_density = pi.density + InitData.tau * pi.right_side_density;
        
        // скорости
//        pi.right_side_velocity = pi.right_side_velocity * 0;
//        for (auto &pj : particles) {
//            r_ij = pi.coordinate - pj.coordinate;
//            if (abs(r_ij) <= InitData.h) {
//                pi.right_side_velocity = pi.right_side_velocity + (pow(pi.density, -1) + pow(pj.density, -1)) * grad_W(pi.coordinate, pj.coordinate, InitData);
//            }
//        }
//        pi.right_side_velocity = pi.right_side_velocity * (-1) * pow(InitData.Cs, 2) * InitData.m_g;
//        pi.right_side_velocity = pi.velocity + InitData.tau * pi.right_side_velocity;
        
        // скорости (прямая аппроксимация градиента давления)
        pi.right_side_velocity = pi.right_side_velocity * 0;
        for (auto &pj : particles) {
            r_ij = pi.coordinate - pj.coordinate;
            if (abs(r_ij) <= InitData.h) {
                pi.right_side_velocity = pi.right_side_velocity + grad_W(pi.coordinate, pj.coordinate, InitData);
            }
        }
        pi.right_side_velocity = pi.right_side_velocity * (-1) * pow(InitData.Cs, 2) * InitData.m_g / pi.density;
        pi.right_side_velocity = pi.velocity + InitData.tau * pi.right_side_velocity;
        
        // давления
//        pi.right_side_pressure = pow(InitData.Cs, 2) * pi.density;
    }
}

void update(std::vector<particle> &particles) {
    // обновление значений искомых величин
    for (auto &pi : particles) {
        pi.coordinate = pi.right_side_coordinate;
        
        pi.velocity = pi.right_side_velocity;
        
        pi.density = pi.right_side_density;
        
//        pi.pressure = pi.right_side_pressure;
    }
}
