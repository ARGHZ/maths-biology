from numpy import *
from matplotlib.pyplot import *
from scipy.integrate import odeint


# divuno=km/(deltam+mu)
# divdos=((L*n)/((Ln)+(Km*n)))
# divtres=(E/(E+ke))
# divcuatro=(L/(L+Kl))

def biestabilidad(s, t, km, kp, deltam, mu, n, Km, deltap, ke, kl, Kl, E):
    P, L = s
    dP = kp * (km / (deltam + mu)) * ((L**n) / ((L**n) + (Km**n))) - (deltap + mu) * P
    dL = (ke * (E / (E + ke)) * P) - (kl * (L / (L + Kl)) * P) - mu * L
    # dP=(kp*divuno*divdos)-(deltap+mu)*P
    # dL=(ke*divtres*P)-(kl*divcuatro*p)-mu*L

    return ([dP, dL])


km = 0.18
kp = 18.1
deltam = 0.46
mu = 0.02
n = 4
Km = 10
deltap = 0.01
ke = 500
kl = 450
Kl = 1400
E = 3
P0 = 0.02
L0 = 0.01
t = linspace(0, 400, 100)

s = odeint(biestabilidad,[P0, L0], t, args=(km, kp, deltam, mu, n, Km, deltap, ke, kl, Kl, E))

P = s[:, 0]
L = s[:, 1]

plot(t, P, label='Promotor')
plot(t, L, label='lactosa')
legend()
show()