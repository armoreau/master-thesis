import numpy as np
import matplotlib.pyplot as plt
from Options import Options
from ode45 import ode45
#20 minutes

def odefcn(x,y) :
    epsilon = 1e-2
    return ((1 - x)*y - y**2)/epsilon

epsilon = 1e-6
y0 = np.array([1])
xspan = np.array([0, 2])


options = Options()
#[x1,y1] = ode15s(@odefcn,xspan,y0,options);

# Impose non-negativity constraint
options.odeset('NonNegative',np.array([0]))
res = ode45(odefcn,xspan,y0,options)

fig = plt.figure()

plt.title('The knee problem');
plt.xlabel('x');
plt.ylabel('solution y')
plt.plot(res.t,res.y[0],label='Non-negativity')
plt.show()
