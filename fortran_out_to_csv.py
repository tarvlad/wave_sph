import sys
import matplotlib.pyplot as plt

if __name__ == '__main__':
    inp = [i.split() for i in open(sys.argv[1], 'r').read().split('\n')]
    out = open(sys.argv[2], 'w')

    tau = float(inp[0][1])
    h = float(inp[1][1])
    n = float(inp[2][1])
    time = float(inp[3][1])

    max_err_v = float(inp[5][1])
    max_err_rho = float(inp[6][1])
    inp = inp[11:]

    data = []
    i = 0
    while i <= len(inp) - 3:
        data.append((
            float(inp[i][0]),
            float(inp[i + 1][0]), float(inp[i + 1][1]),
            float(inp[i + 2][0]), float(inp[i + 2][1]),
        ))
        i += 3

    out.write(f"tau, h, n, time\n")
    out.write(f"{tau}, {h}, {n}, {time}\n")
    out.write("\n")
    out.write("max_err_v, max_err_rho\n")
    out.write(f"{max_err_v}, {max_err_rho}\n")
    out.write("\n")

    out.write("x, v, rho, a_v, a_rho\n")
    for line in data:
        for i in line:
            out.write(str(i) + ",")
        out.write("\n")
    out.close()

    data_v = []
    for i in data:
        data_v.append((i[0], i[1], i[3]))

    data_rho = []
    for i in data:
        data_rho.append((i[0], i[2], i[4]))

    plt.plot(
        [i[0] for i in data_v], [i[1] for i in data_v],
        [i[0] for i in data_v], [i[2] for i in data_v]
    )
    plt.savefig(sys.argv[3])
    plt.clf()

    plt.plot(
        [i[0] for i in data_rho], [i[1] for i in data_rho],
        [i[0] for i in data_rho], [i[2] for i in data_rho]
    )
    plt.savefig(sys.argv[4])