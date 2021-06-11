import matplotlib.pyplot as plt

#Add parent folder to the path. Code taken from https://codeolives.com/2020/01/10/python-reference-module-in-parent-directory/
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from ode45 import ode45

def ode(t,y,beta,gamma,N) :
    return [-beta * y[0] * y[1]/N, beta * y[0] * y[1]/N - gamma * y[1], gamma * y[1]]

tspan = [0, 60]
y0 = [99, 1, 0]
beta = 400/365
gamma = 1/13
N = sum(y0)
argsup = [beta, gamma, N]

sol = ode45(ode, tspan, y0, options = None, varargin = argsup)

#Plot result
fig = plt.figure()
plt.title('Epidemic problem')
plt.xlabel('time [d]')
plt.ylabel('% population')
plt.plot(sol.t,sol.y[0],label="healthy")
plt.plot(sol.t,sol.y[1],label="infected")
plt.plot(sol.t,sol.y[2],label="immune")
plt.legend()
plt.show()