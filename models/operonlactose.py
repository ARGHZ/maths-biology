# coding=utf-8
#MODELO CON VACUNACIÓN
'''
Donde:
    * alfa    --> Nacimientos de la población.
    * beta    --> Probabilidad de contacto entre un suceptible y un infeccioso.
    * gamma   --> Tasa de recuperación (I-->R).
    * delta   --> Tasa natural de muerte independiente del tiempo.
    * epsilon --> Tasa de muertes debido a la enfermedad.
    * R0      --> Número reproductivo básico (Número de contagios debido a un infeccioso mientras es infeccioso)
'''
from numpy import *
from matplotlib.pyplot import *
from mpl_toolkits import mplot3d
from scipy.integrate import odeint


def extracellularlactoseratio(E, K_E):
    result = (E / (E + K_E))

    return result


def intracellularlactoseratio(L, K_L):
    result = (L / (L + K_L))

    return result


def messengersfunction(k_m, gamma_m, mu, L, K_m, n):
    result = (k_m / (gamma_m + mu)) * (L**n / (L **n + K_m**n))

    return result


def proteinsfunction(k_p, gamma_p, mu, k_m, gamma_m, L, K_m, n):
    result = (k_p / (gamma_p + mu)) * (k_m / (gamma_m + mu)) * (L**n / (L**n + K_m**n))

    return result


def extracellularlactosefunctionratio(E, t, k_E, mu, K_E, k_m, k_p, gamma_m, gamma_p):
    numerador = (k_E * E * k_m * k_p)
    denominador = (mu * (E + K_E) * (gamma_m + mu) * (gamma_p + mu))
    result = numerador / denominador

    return result


def lactosaV(L, t, kE, E, kM, kP, n, Mu, KE, KM, GammaM, GammaP):
    numerador = (kE*E*kM*kP*L**n)
    denominador = (Mu * (E + KE) * (GammaM + Mu) * (GammaP + Mu) * ((L**n) + (KM**n)))
    L = numerador / denominador
    return L


# This seems strange
def operonlactoseoriginal(L, t, k_E, mu, E, K_L, k_p, gamma_p, k_m, gamma_m, K_m, n):
    first_division = (k_E / mu)
    second_division = (E / (E + K_E))
    #print '{} * {} / ({} + {}) * ({} + {})'.format(k_m, k_p, gamma_m, mu, gamma_p, mu)
    third_division = ((k_m * k_p) / ((gamma_m + mu) * (gamma_p + mu)))
    fourth_division = ((L**n) / ((L**n) + (K_m**n)))
    #print '{} - {} - {} - {}'.format(first_division, second_division, third_division, fourth_division)
    result = first_division * second_division * third_division * fourth_division
    return (result)


def operonlactosefactorized(L, t, k_E, mu, E, K_L, k_p, gamma_p, k_m, gamma_m, K_m, n):
    numerador = (k_E  * E * k_m * k_p * L**n)
    denominador = (mu * (E + K_E) * (gamma_m + mu) * (gamma_p + mu) * ((L**n) + (K_m**n)))
    result = numerador / denominador

    return result

def intracellularlactosefunction(L, t, k_E, mu, E, K_E, k_L, K_L, k_p, gamma_p, k_m, gamma_m, K_m, n):
    left_section = (k_E / mu) * extracellularlactoseratio(E, K_E) - (k_L / mu) * intracellularlactoseratio(L, K_L)
    right_section = proteinsfunction(k_p=k_p, gamma_p=gamma_p, mu=mu, k_m=k_m, gamma_m=gamma_m, L=L, K_m=K_m, n=n)
    result = left_section * right_section

    return result


def basemodel(init_vector, t, k_m, gamma_m, k_p, gamma_p, E, k_E, K_E, k_L, K_L, mu, K_m, n):
    M, P, L = init_vector
    dM = k_m * ((L**n) / ((L**n) + (K_m**n))) - (gamma_m + mu) * M
    dP = (k_p * M) - (gamma_p + mu) * P
    dL = (k_E * P * (E / (E + K_E))) - (k_L * P * (L / (L + K_L))) - (mu * L)
    return (dM, dP, dL)


def normalbasemodel(init_vector, t, k_m, gamma_m, k_p, gamma_p, E, k_E, K_E, k_L, K_L, mu, K_m, n):
    M, P, L = init_vector
    dM = messengersfunction(k_m, gamma_m, mu, L, K_m, n)
    dP = proteinsfunction(k_p, gamma_p, mu, k_m, gamma_m, L, K_m, n)
    dL = intracellularlactosefunction(L, t, k_E, mu, E, K_E, k_L, K_L, k_p, gamma_p, k_m, gamma_m, K_m, n)
    return (dM, dP, dL)


if __name__ == '__main__':
    t = linspace(0., 1000., 1000)
    threshold = np.repeat(1, t.shape[0])
    init_conditions = (1.0, 1.0, 0.8)

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

    integrated_results = odeint(basemodel, init_conditions, t, args=(k_m, gamma_m, k_p, gamma_p, E, k_E, K_E, k_L, K_L, mu, K_m, n))
    '''integrated_results = odeint(normalbasemodel, init_conditions, t,
               args=(k_m, gamma_m, k_p, gamma_p, E, k_E, K_E, k_L, K_L, mu, K_m, n))
    '''
    M = integrated_results[:, 0]
    P = integrated_results[:, 1]
    L = integrated_results[:, 2]

    ylabel('Concentracion')
    xlabel('Tiempo t')
    title('Modelo Operon Lactosa')
    # Crecimiento de M
    plot(t, M, '-', color="blue", label="Crecimiento de M(t)")
    # Crecimiento de P
    plot(t, P, '-', color="red", label="Crecimiento de P(t)")
    # Crecimiento de L
    plot(t, L, '--', color="green", label="Crecimiento de L(t)")
    legend(loc='center right')
    #xlim((-1,200))
    show()
    clf()

    fig = figure()
    ax = axes(projection='3d')
    # Data for a three-dimensional line
    zline = np.linspace(0, 15, 1000)
    xline = np.sin(zline)
    yline = np.cos(zline)
    ax.set_xlabel('M')
    ax.set_ylabel('P')
    ax.set_zlabel('L')
    ax.scatter3D(M, P, L, c=L, cmap='Blues')
    #show()
    clf()

    # Displaying only Intracellular Lactose concentration
    l = linspace(0., 100., 100)
    threshold = np.repeat(1, t.shape[0])
    init_conditions = (1.0, 1.0, .8)

    results_operon = operonlactoseoriginal(l, 0, k_E, mu, 0.00025, K_L, k_p, gamma_p, k_m, gamma_m, K_m, n)

    title('Modelo Operon Lactosa')
    plot(l, l, '-', color="red", label="L(L)")
    plot(l, results_operon, '-+', color="blue", label="L(t)")
    #xlim((0, 40))
    legend(loc='center right')
    #show()
