# coding=utf-8
from numpy import*
from matplotlib.pyplot import*
from scipy.integrate import odeint

#suma=[]


def vanderpol(r,t,R1,R2,mu1,mu2,mu3):
    x1,y1,x2,y2,x3,y3 = r

    dx1=mu1*(x1-(x1**3)/3-y1)
    dy1=x1/mu1

    dx2=mu2*(x2-(x2**3)/3-y2)+x1/R1
    dy2=x2/mu2+y1/R1

    dx3=mu3*(x3-(x3**3)/3-y3)+x2/R2
    dy3=x3/mu3+y2/R2

    return [dx1, dy1, dx2, dy2, dx3, dy3]


t=linspace(0,100,100)
mu1=3
mu2=5
mu3=7
R1=0.5
R2=15.5
#R1=2
#R2=2



r = odeint(vanderpol,[1,1,1,1,1,1],t,args=(R1,R2,mu1,mu2,mu3,))

#x1=odeint(vanderpol,[[2,2],[1,1],[1,1]],t,args=(R1,R2,mu,))

#x2=odeint(vanderpol,[[2,2],[1,1],[1,1]],t,args=(R1,R2,mu,x1,))

#x3=odeint(vanderpol,[2,2],[1,1],[1,1],t,args=(R1,R2,mu,x2,))

#r=odeint(vanderpol,[2,2],[1,1],[1,1],t,args=(R1,R2,mu1,mu2,mu3,))

#plot(r[:,0],r[:,1])

#plot(r[:,2],r[:,3])

#plot(r[:,4],r[:,5])

#numeros=linspace(0,100,100)


#for i in range(t):
 #   suma=sum(r)
    #plot(numeros,r)



plot(t,r[:,0],label='circuito uno')
plot(t,r[:,2],label='circuito dos')
plot(t,r[:,4],label='circuito tres')
#show()
#plot(r[:,0],r[:,2])
legend()
show()

alpha=1.5
beta=-20
gamma=3

#alpha=0.9
#beta=-0.6
#gamma=0.3

for i in range(len(t)):
    #one=alpha*add(r[:,0],r[:,1])
    two=beta*add(r[:,2],r[:,3])
    three=gamma*add(r[:,4],r[:,5])

    todoJunto=add(two, three)

plot(t,todoJunto)
title("ECG")
show()


plot(t,two)
show()
plot(t,three)
show()
