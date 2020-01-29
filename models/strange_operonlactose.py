# coding=utf-8


from numpy import *
from matplotlib.pyplot import *
from mpl_toolkits import mplot3d
from scipy.integrate import odeint





# This seems strange
def intracellularlactosefunctiondivisions(L, t, k_E, mu, E, K_L, k_p, gamma_p, k_m, gamma_m, K_m, n):
    first_division = (k_E / mu)
    second_division = (E / (E + K_E))
    third_division = ((k_m * k_p) / ((gamma_m + mu) * (gamma_p + mu)))
    fourth_division = ((L**n) / ((L**n) + (K_m**n)))
    result = first_division * second_division * third_division * fourth_division
    return result


def quasistationaryintracellularlactosefunction(L, t, k_E, mu, E, K_L, k_p, gamma_p, k_m, gamma_m, K_m, n):
    numerador = (k_E  * E * k_m * k_p * L**n)
    denominador = (mu * (E + K_E) * (gamma_m + mu) * (gamma_p + mu) * ((L**n) + (K_m**n)))
    result = numerador / denominador

    return result


if __name__ == '__main__':
    t = linspace(0., 1000., 1000)
    threshold = np.repeat(1, t.shape[0])
    init_conditions = (1.0, 1.0, 1.0)

    k_m = 0.18
    gamma_m = 0.46
    k_p = 18.1
    gamma_p = 0.01
    E = 1
    k_E = 500
    K_E = 50
    k_L = 450
    K_L = 1400
    mu = 0.02
    K_m = 10
    n = 4

    # Displaying only Intracellular Lactose concentration
    t = linspace(0., 100., 10000)
    threshold = np.repeat(1, t.shape[0])
    init_conditions = (1.0, 1.0, .8)
    lactose_a = odeint(quasistationaryintracellularlactosefunction, init_conditions[2], t,
                       args=(k_E, mu, E, K_L, k_p, gamma_p, k_m, gamma_m, K_m, n))
    lactose_b = odeint(intracellularlactosefunctiondivisions, init_conditions[2], t,
                       args=(k_E, mu, E, K_L, k_p, gamma_p, k_m, gamma_m, K_m, n))

    ylabel('Distribucion de poblacion (Individuos)')
    xlabel('Tiempo (Minutos)')
    title('Modelo Operon Lactosa')
    plot(t, lactose_b, '-', color="red", label="4 Cocientes L(t)")
    plot(t, lactose_a, '-', color="green", label="1 Cociente L(t)")
    legend(loc='center right')
    show()
