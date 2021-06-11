import numpy as np
import matplotlib.pyplot as plt

#Add parent folder to the path. Code taken from https://codeolives.com/2020/01/10/python-reference-module-in-parent-directory/
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from ode45 import ode45

def ode(t,y) :
    dydt = -2*y + 2*np.cos(t)*np.sin(2*t)
    return dydt

tspan = [0,3]
y0 = [-5,-4,-3,-2,-1,0,1,2,3,4,5]

sol = ode45(ode,tspan,y0)

#Plot ode45 approx
fig = plt.figure()
plt.title('Multiple initial condition')
plt.xlabel('t')
plt.ylabel('y')
for i in range(11):
    plt.plot(sol.t,sol.y[i])
plt.show()