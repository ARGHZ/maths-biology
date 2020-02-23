# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def displaycharts(matrix, t):
    x, y = matrix[:, 0], matrix[:, 1]

    plt.plot(t, x, t, y)
    plt.xlabel('t')
    plt.legend(('x', 'y'))
    plt.show()
    plt.clf()

    # phase portrait
    plt.figure()
    plt.plot(x, y)
    plt.plot(x[0], y[0], 'ro')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    plt.clf()


def periodicdisturbancesimulation(Np, Tp, Hp, mu):
    t = np.array([0])
    x_0 = np.array([0.5, 0])
    x = np.zeros((1, 2))

    x[0] = x_0
    for i in range(Np):
        start, stop, num = i * Tp, (i + 1) * Tp, 100
        t_aux = np.linspace(start, stop, num)
        x_aux = odeint(vanderpol, x_0, t_aux, args=(mu,))
        x_0 = x_aux[-1] + np.array([Hp, 0])

        t = np.concatenate([t, t_aux])
        x = np.concatenate([x, x_aux], axis=0)

    return t, x


def simulatedisturbances():
    mu = 10
    Np = 20
    Hp = 0.75
    Tp = 77

    t, x = periodicdisturbancesimulation(Np, Tp, Hp, mu)

    displaycharts(x, t)


def vanderpol(conditions_vector, t, mu_vector, r_vector):
    dx_vector, dy_vector = [0, 0, 0], [0, 0, 0]

    # Sinus
    x, y = conditions_vector[0], conditions_vector[1]
    dx_vector[0] = mu_vector[0] * (x - (1./3.) * (x**3) - y)
    dy_vector[0] = x / mu_vector[0]

    # Auriculum
    #x, y = conditions_vector[2], conditions_vector[3]
    dx_vector[1] = mu_vector[1] * (x - (1./3.) * (x**3) - y) + x / r_vector[0]
    dy_vector[1] = x / mu_vector[1] + y / r_vector[0]

    # Ventriculum
    #x, y = conditions_vector[4], conditions_vector[5]
    dx_vector[2] = mu_vector[2] * (x - (1./3.) * (x**3) - y) + x / r_vector[1]
    dy_vector[2] = x / mu_vector[2] + y / r_vector[1]

    return dx_vector[0], dy_vector[0], dx_vector[1], dy_vector[1], dx_vector[2], dy_vector[2]


def disturbancesimulation(t, mu, tper, hper, init_conds, r):
    idx = np.where(t <= tper)[0]
    t_0 = t[idx]

    idx = np.where(t > tper)[0]
    t_1 = t[idx]

    #x_0 = odeint(vanderpol, (2, 0), t_0, args=(mu,))
    x_0 = odeint(vanderpol, init_conds, t_0, args=(mu, r,))

    #x_1 = odeint(vanderpol, x_0[-1] + np.array([hper, 0]), t_1, args=(mu,))
    x_1 = odeint(vanderpol, x_0[-1] + np.array([0, 0, hper, 0, 0, 0]), t_1, args=(mu, r,))

    x = np.concatenate([x_0, x_1], axis=0)

    return x


def simulatevanderpol():
    t = np.linspace(0, 1500, 10000)
    params_list = (((100, 100, 100), (0, 0.5, 0, 0.5, 0, 0.5), (0.5, 1)), ((20, 100, 200), (0, 0.5, 0.5, 0, 0.5, 0), (5, 10)), ((30, 150, 300), (0, 0.5, 0.5, 0, 0.5, 0), (5, 10)), ((40, 200, 400), (0, 0.5, 0.5, 0, 0.5, 0), (5, 10)))

    for param_set in params_list[:1]:
        mu, init_conds, r = param_set[0], param_set[1], param_set[-1]

        solutions = disturbancesimulation(t, mu, 479, 1380*3, init_conds, r)
        #solutions = odeint(vanderpol, init_conds, t, args=(mu, r,))

        plt.plot(t, solutions[:, 0], label='Sinus circuit')
        plt.plot(t, solutions[:, 2], label='Auriculum circuit')
        plt.plot(t, solutions[:, 4], label='Ventriculum circuit')
        plt.legend()
        plt.show()
        plt.clf()

        weights = (90, 45, 7)

        '''for ith_t in range(t.shape[0]):
            auriculum = weights[0] * np.add(solutions[:,2], solutions[:,3])
            ventriculum = weights[2] * np.add(solutions[:,4], solutions[:,5])

            sum_up = np.multiply(auriculum, ventriculum)
        '''
        x_vector = solutions[:, 0] * solutions[:, 2] * solutions[:, 4]
        y_vector = solutions[:, 1] * solutions[:, 3] * solutions[:, 5]

        plt.plot(t, x_vector)
        plt.plot(t, y_vector)
        plt.title('ElectroCardioGram')
        plt.show()
        plt.clf()

        plt.plot(x_vector, y_vector)
        plt.title('Phase Space')
        # plt.show()


if __name__ == '__main__':
    simulatevanderpol()
    # simulatedisturbances()
