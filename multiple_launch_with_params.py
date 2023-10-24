import os

if __name__ == '__main__':
    source = open("main.cpp", 'r').read()
    params = (
        (0.8, 63),
    )

    for i in params:
        new_src = source.replace(
            "const int N_gas_particles_per_unit_side = 63;",
            f"const int N_gas_particles_per_unit_side = {i[1]};"
        ).replace(
            "const double h = 0.8;",
            f"const double h = {i[0]};"
        )
        out = open(f"autogen_source/main{i[1]}.cpp")
        out.write(new_src)
        out.close()

    for i in params:
        os.system("clang++ -std=c++20 -Wall -march=native -Ofast -c" +
                  f"autogen_source/main{i[1]}.cpp -o autogen_out/wave_{i[1]}.exe")
        os.system(f"autogen_out/wave_{i[1]}.exe")