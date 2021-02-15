import numpy as np
import matplotlib.pyplot as plt
from ode45 import ode45

def expo(t, y):
    return np.array([-0.5 * y])


########### ODE45 approx

tspan = [10,9,8,5,4,3.33,2.6,1,0]
y0 = [2]
sol = ode45(expo,tspan,y0)

########### Plot ode

fig = plt.figure()
plt.title('Ode45 approx')
plt.xlabel('t')
plt.ylabel('y')
plt.plot(sol.t,sol.y[0])
plt.show()