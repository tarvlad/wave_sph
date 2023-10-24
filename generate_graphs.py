import matplotlib.pyplot as plt
import sys


# Reads data from file in argv[1]
# Out graphs names specified in argv[2], argv[3]
# Input data format:
# x, u, rho, an_u, an_rho
if __name__ == '__main__':
    data = [i.split(',') for i in open(sys.argv[1], 'r').read().rstrip().split('\n')[1:]]

    data_v = [tuple(map(lambda x: float(x), (i[0], i[1], i[3]))) for i in data]
    data_rho = [tuple(map(lambda x: float(x), (i[0], i[2], i[4]))) for i in data]

    plt.plot(
        [i[0] for i in data_v], [i[1] for i in data_v],
        [i[0] for i in data_v], [i[2] for i in data_v]
    )
    plt.savefig(sys.argv[2])
    plt.clf()

    plt.plot(
        [i[0] for i in data_rho], [i[1] for i in data_rho],
        [i[0] for i in data_rho], [i[2] for i in data_rho]
    )
    plt.savefig(sys.argv[3])
    plt.clf()

