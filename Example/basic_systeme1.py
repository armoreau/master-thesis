import matplotlib.pyplot as plt

#Add parent folder to the path. Code taken from https://codeolives.com/2020/01/10/python-reference-module-in-parent-directory/
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from ode45 import ode45

def syst(t, y):
    return [y[1]+y[2], -y[0]+2*y[1]+y[2], y[0]+y[2]]

tspan = [0,2]
y0 = [2,4,6]

sol = ode45(syst,tspan,y0)
    
#Plot ode45 approx
fig = plt.figure()
plt.title('Ode45 approx')
plt.xlabel('t')
plt.ylabel('y')
plt.plot(sol.t,sol.y[0],label="y_1(t)")
plt.plot(sol.t,sol.y[1],label="y_2(t)")
plt.plot(sol.t,sol.y[2],label="y_3(t)")
plt.legend()
plt.show()