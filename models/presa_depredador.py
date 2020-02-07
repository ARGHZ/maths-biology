# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def presadepreadormodel(individuals, t, alpha, beta, gamma, delta, epsilon):
    x, y = individuals

    # alpha: tasa de crecimiento de las presas
    # beta: éxito en la casa del depredador
    # gamma: tasa de decrecimiento de los depredadores
    # delta: éxito en la casa y cuánto alimenta cazar una presa al depredador

    dx = alpha * x - beta * x * y
    dy = - gamma * y + delta * x * y + epsilon * y

    return np.array([dx, dy])


# Parámetros
alpha = 2
beta = 0.8
gamma = 0.4
delta = 0.2
epsilon = 0.5

# Condiciones iniciales
x0 = 10  # Presas
y0 = 2  # Depredadores
conds_iniciales = np.array([x0, y0])

# Condiciones para integración
tf = 10
N = 1000
t = np.linspace(0, tf, N)

solucion = odeint(presadepreadormodel, conds_iniciales, t, args=(alpha, beta, gamma, delta, epsilon))

plt.plot(t, solucion[:, 0], label='Presa')
plt.plot(t, solucion[:, 1], label='Depredadores')
plt.ylabel('Numero de individuos')
plt.xlabel('Tiempo t')
plt.legend()
plt.show()
plt.clf()

x_max = np.max(solucion[:,0]) * 1.05
y_max = np.max(solucion[:,1]) * 1.05
x = np.linspace(0, x_max, 25)
y = np.linspace(0, y_max, 25)
xx, yy = np.meshgrid(x, y)
uu, vv = presadepreadormodel((xx, yy), 0, alpha, beta, gamma, delta, epsilon)
norm = np.sqrt(uu**2 + vv**2)
uu = uu / norm
vv = vv / norm
plt.quiver(xx, yy, uu, vv, norm, cmap=plt.cm.gray)
plt.plot(solucion[:, 0], solucion[:, 1])
plt.ylabel('Depredadores')
plt.xlabel('Presas')
plt.show()
plt.clf()
