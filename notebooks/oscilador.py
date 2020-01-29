# coding=utf-8
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint


def movimientoarmonicosimple(init_conditions, t, k, m, B):
    x, y = init_conditions
    B = 0
    equation = -k * x / m + B / m * y
    return (y, equation)


def movimientoarmonicosimpleamortiguado(init_conditions, t, k, m, B):
    x, y = init_conditions
    equation = -k * x / m + B / m * y
    return (y, equation)


def movimientoarmonicosimpleamortiguadofuerza(init_conditions, t, k, m, B, F_0, w_0):
    x, y = init_conditions
    equation = -k * x / m + B / m * y + F_0 / m * np.sin(w_0 * t)
    return (y, equation)


def displayexperiment(type='simple'):
    k = 4.0  # Constante del Resorte
    m = 1.0  # Masa
    B = -1.0  # Constante de fricción
    F_0 = 2.0
    w_0 = np.sqrt(k / m)

    init_conditions = (0.7, 0.5)
    t = np.arange(0, 10, 0.01)

    if type == 'fuerza':
        result = odeint(movimientoarmonicosimpleamortiguadofuerza, init_conditions, t, args=(k, m, B, F_0, w_0))
    elif type == 'amortiguado':
        result = odeint(movimientoarmonicosimpleamortiguado, init_conditions, t, args=(k, m, B,))
    else:
        result = odeint(movimientoarmonicosimple, init_conditions, t, args=(k, m, B,))

    fig, ax2 = plt.subplots(ncols=1, figsize=(10,10))

    xx, yy = result.T

    plt.xlabel('t')
    plt.plot(t, xx, label="Posicion")
    plt.plot(t, yy, label="Velocidad")
    plt.legend()
    plt.show()

if __name__ == '__main__':
    displayexperiment('fuerza')