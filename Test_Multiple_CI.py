import numpy as np
import matplotlib.pyplot as plt
from ode45 import ode45

def ode(t,y) :
    dydt = -2*y + 2*np.cos(t)*np.sin(2*t)
    return dydt

tspan = [0,3]
y0 = [-3,-2,-1,0,1,2,3]
sol = ode45(ode,tspan,y0)

fig = plt.figure()
plt.title('Multiple condition initial')
plt.xlabel('t')
plt.ylabel('y')
for i in range(7):
    plt.plot(sol.t,sol.y[i])
plt.show()
