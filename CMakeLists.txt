cmake_minimum_required(VERSION 3.5)

enable_language(Fortran)
add_compile_options("-fdefault-real-8")
add_compile_options("-std=f2018")


set(SOURCE wave_sph.f90)
add_executable(wave_1d ${SOURCE})