import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def exponential(x, t, gamma_param):
    dx = gamma_param * x
    return dx


def competition(r, t, a, b, rho):
    x, y = r
    dx = x * (1 - x - a * y)
    dy = rho * y * (1 - y - b * x)

    return (dx, dy)


def logispopulation():
    pass


if __name__ == '__main__':
    gamma_param, x_0 = 1, 1
    t = np.linspace(0, 5, 10)
    x = odeint(exponential, x_0, t, args=(gamma_param, ))
    '''plt.plot(t, x)
    plt.xlabel('Time t')
    plt.ylabel('x values')
    plt.plot(t, x_0 * np.exp(gamma_param * t))
    plt.show()'''

    a, b= 2.0, 2.0
    rho = 4
    t = np.linspace(0.0, 9.0, 100.0)
    r = odeint(competition, (1, 1), t, args=(a, b, rho))
    x = r[:, 0]
    y = r[:, 1]
    #plt.plot(t, x, '*')
    plt.xlabel('Time t')
    plt.plot(x, y, 'o')
    r = odeint(competition, (1, 0.9), t, args=(a, b, rho))
    x = r[:, 0]
    y = r[:, 1]
    plt.plot(x, y, 'o')
    '''r = odeint(competition, (0.85, 1), t, args=(a, b, rho))
    x = r[:, 0]
    y = r[:, 1]
    plt.plot(x, y)
    r = odeint(competition, (0.6, 1), t, args=(a, b, rho))
    x = r[:, 0]
    y = r[:, 1]
    plt.plot(x, y)
    r = odeint(competition, (0.4, 1), t, args=(a, b, rho))
    x = r[:, 0]
    y = r[:, 1]
    plt.plot(x, y)'''
    plt.show()

