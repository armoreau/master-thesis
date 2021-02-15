import ode45
import Options
import numpy as np
import matplotlib.pyplot as plt

def odefcn(t,y,A,B) :
    return np.array([y[1], (A/B)*t* y[0]])

######### Modify options

tspan = np.array([0, 5])
y0 = np.array([0, 0.01])
A=1
B=2

myOptions = Options.Options() #Create an option with default value.
myOptions.odeset('RelTol',np.array([1e-5]))
myOptions.odeset('AbsTol',1e-8*np.ones(y0.size))
myOptions.odeset('Refine',1)
myOptions.odeset('NormControl',False)
myOptions.odeset('MaxStep',np.array([10]))
myOptions.odeset('InitialStep',np.array([0.1]))

sol = ode45.ode45(odefcn,tspan,y0, myOptions, varargin = (A,B))

########### Plot ode

fig = plt.figure()
plt.title('Ode45 approx')
plt.xlabel('t')
plt.ylabel('y')
plt.plot(sol.t,sol.y[0],label="y_1(t)")
plt.plot(sol.t,sol.y[1],label="y_2(t)")
plt.legend()
plt.show()