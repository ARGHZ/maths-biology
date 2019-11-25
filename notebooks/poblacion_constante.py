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
from scipy.integrate import odeint


def basemodel(r_vector, t, alfa, beta, gamma, delta, epsilon, V):
    S, I, R = r_vector
    dS = -beta * S * I
    dI = (beta * S * I) - (I * gamma)
    dR = gamma * I
    return (dS, dI, dR)


def normalbasemodel(r_vector, t, r_0):
    s_script, i_script = r_vector
    dS = -r_0 * s_script * i_script
    dI = (r_0 * s_script * i_script) - i_script
    return (dS, dI)


def vaccinebasemodel(r_vector, t, beta, gamma, eta):
    S, I, R = r_vector
    dS = (-beta * S * I) - (eta * S)
    dI = (beta * S * I) - (I * gamma)
    dR = (gamma * I) + (eta * S)
    return (dS, dI, dR)


def normalvaccinebasemodel(r_vector, t, r_0, tau):
    s_script, i_script = r_vector
    dS = (-r_0 * i_script - tau) * s_script
    dI = (r_0 * s_script - 1) - i_script
    return (dS, dI)


def experimentbasemodels():
    t = linspace(0., 2.5, 10000)
    threshold = np.repeat(1, t.shape[0])
    init_conditions = (500000, 500000, 0.0)

    alfa = 12500
    beta = 2 / 8.e7
    gamma = 10.
    delta = .0125
    epsilon = 0
    eta = 0.1
    V = 0.0

    E = odeint(basemodel, init_conditions, t, args=(alfa, beta, gamma, delta, epsilon, V))
    S = E[:, 0]
    I = E[:, 1]
    R = E[:, 2]

    ylabel('Distribucion de poblacion (Individuos)')
    xlabel('Tiempo (Agnios)')
    title('Modelo SIR con vacunacion\n(alfa={}, beta={}, gamma={},\n delta={}, epsilon={}, V={})'.format(alfa, beta,
                                                                                                         gamma, delta,
                                                                                                         epsilon, V))
    # Crecimiento de S
    plot(t, S, '-', color="blue", label="Crecimeinto de S(t)")
    # Crecimiento de I
    plot(t, I, '-', color="red", label="Crecimeinto de I(t)")
    # Crecimiento de R
    plot(t, R, '-', color="green", label="Crecimeinto de R(t)")
    plot(t, threshold, 'm--')
    legend(loc='center right')
    show()

    E = odeint(vaccinebasemodel, init_conditions, t, args=(beta, gamma, eta))
    S = E[:, 0]
    I = E[:, 1]
    R = E[:, 2]

    ylabel('Distribucion de poblacion (Individuos)')
    xlabel('Tiempo (Agnios)')
    title('Modelo SIR con vacunacion\n(alfa={}, beta={}, gamma={},\n delta={}, epsilon={}, V={})'.format(alfa, beta,
                                                                                                         gamma, delta,
                                                                                                         epsilon, V))
    # Crecimiento de S
    plot(t, S, '-', color="blue", label="Crecimeinto de S(t)")
    # Crecimiento de I
    plot(t, I, '-', color="red", label="Crecimeinto de I(t)")
    # Crecimiento de R
    plot(t, R, '-', color="green", label="Crecimeinto de R(t)")
    plot(t, threshold, 'm--')
    legend(loc='center right')
    show()


if __name__ == '__main__':
    experimentbasemodels()
