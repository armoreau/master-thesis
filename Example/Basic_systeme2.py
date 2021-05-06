import numpy as np
import matplotlib.pyplot as plt

#Add parent folder to the path. Code taken from https://codeolives.com/2020/01/10/python-reference-module-in-parent-directory/
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from ode45 import ode45

#Van der pol equation
def syst(t,y) :
    dydt = np.zeros(2)
    dydt[0] = y[1]
    dydt[1] = (1-y[0]**2)*y[1]-y[0]
    return dydt

########### ODE45 approx

tspan = [0,20]
y0 = [2,0]
sol = ode45(syst,tspan,y0)

########### Plot ode

fig = plt.figure()
plt.title('Ode45 approx')
plt.xlabel('t')
plt.ylabel('y')
plt.plot(sol.t,sol.y[0],label="y_1(t)")
plt.plot(sol.t,sol.y[1],label="y_2(t)")
plt.legend()
plt.show()

