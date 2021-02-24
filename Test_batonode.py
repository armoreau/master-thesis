import numpy as np
import matplotlib.pyplot as plt
from Options import Options
from ode45 import ode45

m1 = 0.1
m2 = 0.1
L = 1
g = 9.81

def dydt(t,y) :   
    dydt = np.array([y[1], m2*L*(y[5]**2)*np.cos(y[4]),y[3],m2*L*(y[5]**2)*np.sin(y[4])-(m1+m2)*g,y[5],-g*L*np.cos(y[4])])
    return dydt

def mass(t,y) :   
    M = np.zeros([6,6])
    M[0,0] = 1
    M[1,1] = m1 + m2
    M[1,5] = -m2*L*np.sin(y[4])
    M[2,2] = 1
    M[3,3] = m1 + m2
    M[3,5] = m2*L*np.cos(y[4])
    M[4,4] = 1
    M[5,1] = -L*np.sin(y[4])
    M[5,3] = L*np.cos(y[4])
    M[5,5] = L**2
    return M
    

tspan = np.linspace(0,4,25)
y0 = [0, 4, 2, 20, -np.pi/2, 2]
options = Options()
options.odeset('Mass',mass)
res = ode45(dydt,tspan,y0,options)


fig = plt.figure()
plt.title('A thrown baton problem with mass matrix M(t,y), solved by ODE45')
plt.xlabel('x')
plt.ylabel('y')

for j in range(len(res.t)) :
    theta = res.y[4,j]
    X = res.y[0,j]
    Y = res.y[2,j]
    xvals = np.array([X, X+L*np.cos(theta)])
    yvals = np.array([Y, Y+L*np.sin(theta)])

    plt.plot(xvals,yvals,xvals[0],yvals[0],'ro',xvals[1],yvals[1],'go')

plt.show()