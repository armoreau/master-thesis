import numpy as np
import matplotlib.pyplot as plt

#Add parent folder to the path. Code taken from https://codeolives.com/2020/01/10/python-reference-module-in-parent-directory/
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from ode45 import ode45
from odeoptions import Odeoptions

def sin(t,y) :
    dydt = np.cos(t)
    return dydt

tspan = [0,10]
y0 = [0]

myOptions = Odeoptions() #Create an option with default value.
myOptions.odeset('RelTol',1e-5)
myOptions.odeset('AbsTol',1e-8)
myOptions.odeset('Refine',10)
myOptions.odeset('NormControl',True)
myOptions.odeset('MaxStep',1)
myOptions.odeset('InitialStep',0.1)

sol = ode45(sin,tspan,y0, myOptions)

#Plot ode45 approx
fig = plt.figure()
plt.title('Ode45 approx')
plt.xlabel('t')
plt.ylabel('y')
plt.plot(sol.t,sol.y[0])
plt.show()