import numpy as np
import matplotlib.pyplot as plt
from Options import Options
from ode45 import ode45

def dydt(t,y) :
    # Derivative function -- mu and mustar shared with the outer function.
    mu = 1 / 82.45
    mustar = 1 - mu
      
    r13 = ((y[0] + mu)**2 + y[1]**2) ** 1.5
    r23 = ((y[0] - mustar)**2 + y[1]**2) ** 1.5
    return np.array([ y[2],y[3],2*y[3] + y[0] - mustar*((y[0]+mu)/r13) - mu*((y[0]-mustar)/r23),-2*y[2] + y[1] - mustar*(y[1]/r13) - mu*(y[1]/r23)])
    
def events(t,y):
    y0 = np.array([1.2, 0, 0, -1.04935750983031990726])
    dDSQdt = 2 * np.dot((y[0:2]-y0[0:2]),y[2:4])
    value = np.array([dDSQdt, dDSQdt])
    isterminal = np.array([1, 0]) # stop at local minimum
    direction = np.array([1, -1]) # [local minimum, local maximum]
    return [value,isterminal,direction]


y0 = np.array([1.2, 0, 0, -1.04935750983031990726])
tspan = np.array([0.0, 7.0])
options = Options()
options.odeset('Events',events)
options.odeset('RelTol',np.array([1e-5]))
options.odeset('AbsTol',np.ones([4])*1e-4)

res = ode45(dydt,tspan,y0,options)

fig = plt.figure()
plt.title('Restricted three body problem')
plt.xlabel('x(t)')
plt.ylabel('y(t)')
plt.plot(res.y[0],res.y[1])
plt.show()