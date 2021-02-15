import numpy as np
import matplotlib.pyplot as plt
from ode45 import ode45

def syst(t,y) :
    dydt = np.zeros(2)
    dydt[0] = y[1]
    dydt[1] = (5-y[1]**2)*y[1]-y[0]
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

