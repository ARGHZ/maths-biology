import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def susceptiblediffunction(s_individuals, t, alpha, beta, delta, upsilon, i_individuals):
    ds = (1 - upsilon) * alpha - beta * (s_individuals * i_individuals) - delta * s_individuals
    return ds


def infectioningdiffunction(i_individuals, t, gamma, beta, delta, epsilon, s_individuals):
    ds = beta * s_individuals * i_individuals - (gamma + delta + epsilon) * i_individuals
    return ds


def resistancediffunction(r_individuals, t, gamma, delta, upsilon, alpha, i_individuals):
    ds = gamma * i_individuals - delta * r_individuals + upsilon * alpha
    return ds


def simodel(init_conditions, t, alpha, beta, gamma, delta, epsilon, upsilon):
    s, i, r = init_conditions
    ds = susceptiblediffunction(s, t, alpha, beta, delta, upsilon, i)
    di = infectioningdiffunction(i, t, gamma, beta, delta, epsilon, s)
    dr = resistancediffunction(r, t, gamma, delta, upsilon, alpha, i)

    return (ds, di, dr)


def plotsir(t, alpha = 0.1, beta = 0.1, gamma = 0.1, delta = 0.2, epsilon = 0.2, upsilon = 0.01):
    init_s = np.arange(0.1, 1.01, 0.1)
    init_i = np.arange(0.1, 1.01, 0.1)

    init_population = []
    for step_s in init_s:
        for step_y in init_i:
            pair_values = (step_s, step_y)
            init_population.append(pair_values)

    #for s_0, i_0 in init_population:
    init_conditions = (10, 0, 0)
    process_result = odeint(simodel, init_conditions, t, args=(alpha, beta, gamma, delta, epsilon, upsilon))
    s = process_result[:, 0]
    i = process_result[:, 1]
    r = process_result[:, 1]
    plt.plot(s, 'g', i, 'r', r, 'c')
    #plt.title("S: {} | I: {} | R: {}".format(init_conditions[0], init_conditions[1], init_conditions[2]))
    #plt.xlabel('Population Y')
    #plt.ylabel('Population X')
    plt.show()
    plt.clf()


def experimentsirmodel():
    alpha = 0.1
    beta = 0.1
    gamma = 0.1
    delta = 0.2
    epsilon = 0.2
    upsilon = 0.01

    t = np.linspace(0.0, 365, 1.0)
    plotsir(t, alpha, beta, gamma, delta, epsilon, upsilon)


if __name__ == '__main__':
    # plotexponential()
    # plotlogistic()
    # experimentlotkavolterra()
    experimentsirmodel()
