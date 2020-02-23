from numpy import *
from scipy.integrate import odeint
import matplotlib.pylab as plt


def vanderpol(r, t, mu):
	x, y = r
	dx = mu * (x - (x**3)/3 -y)
	dy = x/mu
	return array([dx, dy])


def perturbation(t, mu, tper, hper):
	idx = where(t <= tper)[0]
	t0 = t[idx]
	idx = where(t > tper)[0]
	t1 = t[idx]
	x0 = odeint(vanderpol, [2, 0], t0, args=(mu,))
	x1 = odeint(vanderpol, x0[-1] + array([hper, 0]), t1, args=(mu,))
	x = concatenate([x0, x1], axis=0)
	return x


def periodicPerturbation(Np, Tp, Hp, mu):
	t = array([0])
	x0 = array([0, 2])
	x = zeros((1, 2))
	x[0] = x0
	for i in range(Np):
		taux = linspace(i*Tp, (i+1)*Tp, 100)
		xaux = odeint(vanderpol, x0, taux, args=(mu,))
		x0 = xaux[-1] + array([Hp, 0])
		t = concatenate([t, taux])
		x = concatenate([x, xaux], axis=0)
	return t, x


if __name__ == '__main__':
	mu = 10
	Np = 20
	Hp = 0.75
	Tp = 17
	t, x = periodicPerturbation(Np, Tp, Hp, mu)
	plt.plot(x[:,0], x[:,1])
	plt.show()
