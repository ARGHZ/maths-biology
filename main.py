from numpy import*
from matplotlib.pyplot import*
from scipy.integrate import odeint

def depredadorpresa(r,t,alpha,beta,gamma,delta,epsilon):
    x,y=r
    dx=alpha*x-beta*x*y
    dy=epsilon*y+gamma*x*y-delta*y
    return (dx,dy)

def plotdepredadorpresa(alpha, beta, gamma, delta, epsilon, t):

    r=odeint(depredadorpresa,[3,1],t,args=(alpha,beta,gamma,delta,epsilon))

    x=r[:,0]
    y=r[:,1]

    plot(t,x,label='Presa')
    plot(t,y,label='Depredador')

    legend()
    show()


def plotphase(alpha, beta, gamma, delta, epsilon, t):
    linspace(0,100,100000)

    init_conditions = ((3, 2), (3, 3), (3, 1), (5,2), (4, 1))

    for init_values in init_conditions:
        r=odeint(depredadorpresa, init_values, t, args=(alpha,beta,gamma,delta,epsilon))
        x=r[:,0]
        y=r[:,1]

        plot(x, alpha-beta*y)

    legend()
    show()


if __name__ == '__main__':
    alpha, beta, gamma, delta, epsilon = 2, 1, 0.4, 1, 0.01
    t = linspace(0, 10, 100)

    plotdepredadorpresa(alpha, beta, gamma, delta, epsilon, t)
    plotphase(alpha, beta, gamma, delta, epsilon, t)

