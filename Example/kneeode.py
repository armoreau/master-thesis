import numpy as np
import matplotlib.pyplot as plt

#Add parent folder to the path. Code taken from https://codeolives.com/2020/01/10/python-reference-module-in-parent-directory/
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from odeoptions import Odeoptions
from ode45 import ode45

def odefcn(x,y) :
    epsilon = 1e-2
    return ((1 - x)*y - y**2)/epsilon

epsilon = 1e-6
y0 = [1]
xspan = [0, 2]
options = Odeoptions()
options.odeset('NonNegative',[0])

res = ode45(odefcn,xspan,y0,options)

fig = plt.figure()

plt.title('The knee problem');
plt.xlabel('x');
plt.ylabel('y')
plt.plot(res.t,res.y[0],label='Non-negativity')
plt.show()
