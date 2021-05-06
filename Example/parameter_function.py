import numpy as np
import matplotlib.pyplot as plt

#Add parent folder to the path. Code taken from https://codeolives.com/2020/01/10/python-reference-module-in-parent-directory/
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from ode45 import ode45

def odefcn(t,y,A,B) :
    return np.array([y[1], (A/B)*t* y[0]])

tspan = [0, 5]
y0 = [0, 0.01]
A=1
B=2

sol = ode45(odefcn, tspan, y0, options = None, varargin = (A,B))
fig = plt.figure()
plt.title('parameter function')
plt.xlabel('t')
plt.ylabel('y')
plt.plot(sol.t,sol.y[0],label="y_1")
plt.plot(sol.t,sol.y[1],label="y_2")
plt.legend()
plt.show()