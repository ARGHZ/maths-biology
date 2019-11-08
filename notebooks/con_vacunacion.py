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


def SIRV(E,t,alfa,beta,gamma,delta,epsilon,V):
    S,I,R=E
    dS=(1-V)*alfa-beta*S*I-delta*S
    dI=beta*S*I-(gamma+delta+epsilon)*I
    dR=gamma*I-delta*R+V*alfa
    return[dS,dI,dR]

t = linspace(0.,2.5,10000)
threshold = np.repeat(1,t.shape[0])

init_conditions = (500000, 500000, 0.0)

alfa=12500
beta=2/8.e7
gamma=10.
delta=.0125
epsilon=0
V=0.0

E=odeint(SIRV,init_conditions,t,args=(alfa,beta,gamma,delta,epsilon,V))
S=E[:,0]
I=E[:,1]
R=E[:,2]

ylabel('Distribucion de poblacion (Individuos)')
xlabel('Tiempo (Agnios)')
title('Modelo SIR con vacunacion\n(alfa={}, beta={}, gamma={},\n delta={}, epsilon={}, V={})'.format(alfa,beta,gamma,delta,epsilon,V))

#Crecimiento de S
plot(t,S, '-', color="blue", label="Crecimeinto de S(t)")
#Crecimiento de I
plot(t,I, '-', color="red", label="Crecimeinto de I(t)")
#Crecimiento de R
plot(t,R, '-', color="green", label="Crecimeinto de R(t)")

plot(t, threshold, 'm--')

legend(loc='center right')

show()