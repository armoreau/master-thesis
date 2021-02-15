import numpy as np
import matplotlib.pyplot as plt
from ode45 import ode45

def sin(t,y) :
    dydt = np.cos(t)
    return dydt

def sol_sin(t) : #condition init y[0] = 1
    return np.sin(t)

########### ODE45 approx

tspan = [0,10]
y0 = [0]
sol = ode45(sin,tspan,y0)

########### True sol

true_sol = np.zeros([1,len(sol.t)])
for i in range(len(sol.t)) :
    true_sol[:,i] = sol_sin(sol.t[i])

########### Plot ode
    
fig = plt.figure()
plt.title('Ode45 approx')
plt.xlabel('t')
plt.ylabel('y')
plt.plot(sol.t,sol.y[0])

########### Plot True sol
fig = plt.figure()
plt.title('True sol')
plt.xlabel('t')
plt.ylabel('y')
plt.plot(sol.t,true_sol[0])
plt.show()
