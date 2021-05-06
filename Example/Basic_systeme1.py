import numpy as np
import matplotlib.pyplot as plt

#Add parent folder to the path. Code taken from https://codeolives.com/2020/01/10/python-reference-module-in-parent-directory/
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from ode45 import ode45

def syst(t, y):
    return np.array([y[1]+y[2], -y[0]+2*y[1]+y[2], y[0]+y[2]])
def sol_syst(t) : #CI : y[0](0) = 2, y[1](0) = 4, y[2](0) = 6,
    return np.array([-3+5*np.exp(2*t), -3+2*np.exp(t)+5*np.exp(2*t), 3-2*np.exp(t)+5*np.exp(2*t)])

########### ODE45 approx

tspan = [0,2]
y0 = [2,4,6]
sol = ode45(syst,tspan,y0)

########### True sol

true_sol = np.zeros([3,len(sol.t)])
for i in range(len(sol.t)) :
    true_sol[:,i] = sol_syst(sol.t[i])
    
########### Plot ode

fig = plt.figure()
plt.title('Ode45 approx')
plt.xlabel('t')
plt.ylabel('y')
plt.plot(sol.t,sol.y[0],label="y_1(t)")
plt.plot(sol.t,sol.y[1],label="y_2(t)")
plt.plot(sol.t,sol.y[2],label="y_3(t)")
plt.legend()

########### Plot True sol

fig = plt.figure()
plt.title('True sol')
plt.xlabel('t')
plt.ylabel('y')
plt.plot(sol.t,true_sol[0],label="y_1(t)")
plt.plot(sol.t,true_sol[1],label="y_2(t)")
plt.plot(sol.t,true_sol[2],label="y_3(t)")
plt.legend()
plt.show()
