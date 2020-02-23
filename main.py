from scipy import linspace
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def vdp(t, z):
    x, y = z
    return [y, mu * (1 - x ** 2) * y - x]


a, b = 0, 180

mus = (10, 50, 100)
styles = ["-", "--", ":"]
t = linspace(a, b, 500)

for mu in mus:
    sol = solve_ivp(vdp, [a, b], [1, 0], t_eval=t)
    plt.plot(sol.y[0], sol.y[1])
    # plt.plot(sol.t, sol.y[0])

# make a little extra horizontal room for legend
plt.legend(["mu={}".format(m) for m in mus])
plt.show()
plt.clf()
