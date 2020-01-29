# coding=utf-8
'''
Modelo Operón Lactosa
M -->
P -->
L --> Lactosa Intracelular
E --> Lactosa Extracelular

kM --> # 1/min
kp --> # 1/min
kL=450 # 1/min
kE=500 # 1/min

KM=10 # micro M
KL=1400 # micro M
KE=50 # micro M

GammaM=0.46 # 1/min
Gammap=0.01 # 1/min
Mu=0.02 # 1/min

n=4 # (Adimensional)

'''

from numpy import *
from matplotlib.pyplot import *
from scipy.integrate import odeint

from mpl_toolkits import mplot3d


#MODELO SIR (sin vacunación)

def operonlactosa(opelac, t, kM, kp, kL, kE, KM, KL, KE, GammaM, Gammap, Mu, E, n):
    M,P,L=opelac
    dM=((kM*(L**n))/((L**n)+(KM*n)))-(GammaM+Mu)*M
    dP=(kp*M)-((Gammap+Mu)*P)
    dL=((kE*E*P)/(E+KE))-(kL*L*(P-Mu*L)/(L+KL))
    return[dM,dP,dL]

M0=1   #
P0=1   #
L0=1   #Lactosa Intracelular
E=1   #Lactosa Extracelular

t=linspace(0.,10,10000)

kM=0.18 # 1/min
kp=18.1 # 1/min
kL=450 # 1/min
kE=500 # 1/min

KM=10 # micro M
KL=1400 # micro M
KE=50 # micro M

GammaM=0.46 # 1/min
Gammap=0.01 # 1/min
Mu=0.02 # 1/min

n=4 # (Adimensional)


opelac=odeint(operonlactosa,[M0,P0,L0],t,args=(kM, kp, kL, kE, KM, KL, KE, GammaM, Gammap, Mu, E, n))
M=opelac[:,0]
P=opelac[:,1]
L=opelac[:,2]

ylabel('Concentracion(Micro-Molar)')
xlabel('Tiempo(Minutos)')
title('Modelo Operon Lactosa\n(M0={}, P0={}, L0={}, E={}, \nkM={}, kp={}, kL={}, kE={}, \nKM={}, KL={}, KE={}, GammaM={}, \nGammap={}, Mu={}, E={}, n={})'.format(M0, P0, L0, E, kM, kp, kL, kE, KM, KL, KE, GammaM, Gammap, Mu, E, n))

#Aumento de Concentración de M
plot(t,M, '-', color="blue", label="Aumento molar de M(t)")
#Aumento de Concentración de P
plot(t,P, '-', color="red", label="Aumento molar de P(t)")
#Aumento de Concentración de L
plot(t,L, '-', color="green", label="Aumento molar de L(t)")

grid()
legend(loc='center right')

show()


ax = axes(projection='3d')
ax.scatter3D(M, P, L, c=L);
ax.set_title('Diagrama de fase');
show()