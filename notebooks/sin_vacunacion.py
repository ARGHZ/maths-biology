# coding=utf-8
#MODELO SIN VACUNACIÓN
# Suceptibles, Infectados, Recuperados (SIR).
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


def SIR(e,t,alfa,beta,gamma,delta,epsilon):
    S,I,R=e
    dS=alfa+beta*S*I+delta*S
    dI=beta*S*I-(gamma+delta+epsilon)*I
    dR=gamma*I-delta*R
    return[dS,dI,dR]

t=linspace(0.,700.,100)

S0=5000
I0=5000
R0=0

alfa=1
beta=0.001
gamma=50
delta=.0125
epsilon=0


e=odeint(SIR,[S0,I0,R0],t,args=(alfa,beta,gamma,delta,epsilon))
S=e[:,0]
I=e[:,1]
R=e[:,2]

ylabel('Distribucion de poblacion (Individuos)')
xlabel('Tiempo (Agnios)')
title('Modelo SIR sin vacunacion\n(alfa={}, beta={}, gamma={},\n delta={}, epsilon={})'.format(alfa,beta,delta,gamma,delta,epsilon))

#Crecimiento de S
plot(t,S, '-', color="blue", label="Crecimeinto de S(t)")
#Crecimiento de I
plot(t,I, '-', color="red", label="Crecimeinto de I(t)")
#Crecimiento de R
plot(t,R, '-', color="green", label="Crecimeinto de R(t)")

legend(loc='center right')

show()