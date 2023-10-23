import sys
import os


def change_time_in_fortran_code(time, i):
    data = open("wave_sph.f90", 'r').read()
    data = data.replace(
        "real, parameter :: time_moment = 1.0",
        "real, parameter :: time_moment = " + str(time) + "!__TIME_TO_INJECT__"
    )
    open(f"autogen_code/wave_sph_{i}.f90", 'w').write(data)


if __name__ == '__main__':
    times = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for i in range(len(times)):
        change_time_in_fortran_code(times[i], i)
        os.system("flang -std=f2018 -fdefault-real-8 " + f"autogen_code/wave_sph_{i}.f90" + " -o wave_1d.exe")
        os.system("wave_1d.exe > out.txt")
        os.system(f"python fortran_out_to_csv.py out.txt out.csv autogen_out/v_{i}.png autogen_out/rho_{i}.png")