import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def exponential(x, t, gamma_param):
    dx = gamma_param * x
    return dx


def nondifexponential(n_0, alpha, t):
    result = n_0 * np.exp(alpha * t)

    return result


def plotexponential():
    gamma_param, x_0 = 1, 1
    t = np.linspace(0, 5, 10)
    x = odeint(exponential, x_0, t, args=(gamma_param, ))
    plt.plot(t, x, '-')
    plt.xlabel('Time t')
    plt.ylabel('x values')
    plt.plot(t, nondifexponential(x_0, gamma_param, t), 'r+')
    plt.show()


def competition(r, t, a, b, rho):
    x, y = r
    dx = x * (1 - x - a * y)
    dy = rho * y * (1 - y - b * x)

    return (dx, dy)


def plotloktavolterra(aggresiveness_x_y, aggresiveness_y_x, rho):
    a, b = aggresiveness_x_y, aggresiveness_y_x

    t = np.linspace(0.0, 9.0, 100.0)
    init_x = np.arange(0.1, 1.01, 0.1)
    init_y = np.arange(0.1, 1.01, 0.1)

    init_population = []
    for x_step in init_x:
        for y_step in init_y:
            pair_values = (x_step, y_step)
            init_population.append(pair_values)

    for x_0, y_0 in init_population:
        r = odeint(competition, (x_0, y_0), t, args=(a, b, rho))
        x = r[:, 0]
        y = r[:, 1]
        plt.plot(x, y, '--')

    plt.title("a: {} | b: {} | rho: {}".format(a, b, rho))
    plt.xlabel('Population Y')
    plt.ylabel('Population X')
    plt.show()
    plt.clf()


def experimentlotkavolterra():

    a_vector = np.arange(1, 4, 1)
    b_vector = np.arange(1, 4, 1)
    rho_vector = np.arange(1, 31, 20)

    for a in a_vector:
        for b in b_vector:
            for rho in rho_vector:
                plotloktavolterra(a, b, rho)


def logistic(n, t, alpha, k):
    alpha = alpha * 1.0
    k = k * 1.0
    result = alpha * (1 - (n/ k)) * n
    return  result


def nondiflogistic(n_0, t, alpha_0, k):
    n_0 = n_0 * 1.0
    alpha_0 = alpha_0 * 1.0
    k = k * 1.0
    numerator = 1
    denominator = (1 / k) + ((1 / n_0) - (1 / k)) * (np.exp(-alpha_0 * t))
    result = numerator / denominator

    return result


def plotlogistic():
    x_0, alpha, k = 0.3, 2, 5
    t = np.linspace(0, 5, 10)
    x = odeint(logistic, x_0, t, args=(alpha, k))
    plt.plot(t, x, '-')
    plt.xlabel('Time t')
    plt.ylabel('x values')
    plt.plot(t, nondiflogistic(x_0, t, alpha, k), 'r+')
    plt.show()


def susceptiblediffunction(r_0, t, script_s, script_i):
    ds = 1 - (r_0 * script_s * script_i) - script_s
    return ds


def infectioningdiffunction(r_0, t, script_s, script_i, rho):
    ds = rho * (r_0 * script_s * script_i - script_i)
    return ds


def simodel(scripts_vector, t, r_0, rho):
    s, i = scripts_vector
    ds = susceptiblediffunction(r_0, t, s, i)
    di = infectioningdiffunction(r_0, t, s, i, rho)

    return (ds, di)


def plotsir(ar_zero, rho):
    r_0, rho = ar_zero, rho
    t = np.linspace(0.0, 9.0, 100.0)
    init_s = np.arange(0.1, 1.01, 0.1)
    init_i = np.arange(0.1, 1.01, 0.1)

    init_population = []
    for step_s in init_s:
        for step_y in init_i:
            pair_values = (step_s, step_y)
            init_population.append(pair_values)

    for s_0, i_0 in init_population:
        s = odeint(simodel, (s_0, i_0), t, args=(r_0, rho))
        x = s[:, 0]
        y = s[:, 1]
        plt.plot(x, y, '--')

    plt.title("a: {} | b: {} | rho: {}".format(a, b, rho))
    plt.xlabel('Population Y')
    plt.ylabel('Population X')
    plt.show()
    plt.clf()


def experimentsirmodel():
    vector_s = np.arange(0, 1.01, 0.5)
    vector_i = np.arange(0, 1.01, 0.5)
    vector_r_0 = np.arange(0, 1.01, 0.5)
    vector_rho = np.arange(0, 1.01, 0.3)

    t = np.linspace(0.0, 9.0, 100.0)

    for s in vector_s:
        for i in vector_i:
            for r_0 in vector_r_0:
                for rho in vector_rho:
                    print "s: {} | i: {} | R0: {} | rho: {}".format(s, i, r_0, rho)



if __name__ == '__main__':
    # plotexponential()
    # plotlogistic()
    # experimentlotkavolterra()
    experimentsirmodel()
