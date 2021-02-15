from scipy.integrate import solve_ivp
import Compare_data
import ode45
import numpy as np
import matplotlib.pyplot as plt
import time


def eq1(t, y):
    return np.array([y[1]+y[2], -y[0]+2*y[1]+y[2], y[0]+y[2]])
def sol1(t) : #CI : y[0](0) = 2, y[1](0) = 4, y[2](0) = 6,
    return np.array([-3+5*np.exp(2*t), -3+2*np.exp(t)+5*np.exp(2*t), 3-2*np.exp(t)+5*np.exp(2*t)])

def eq2(t, y):
    return [-0.5 * y]
def sol2(t) : #condition init y[0] = 2
    return np.array([2*np.exp((-1/2)*t)])

def eq3(t, y):
    return -np.cos(1/y[0])/(y**2)
def sol3(t) : #condition init y[1] = 0
    return np.sin(1/t) - np.sin(1)

def odefcn(t,y,A,B) :
    return np.array([y[1], (A/B)*t* y[0]])

def ode1(t,y) :
    dydt = np.cos(t)
    return dydt

def ode2(t,y) :
    dydt = np.zeros(3)
    dydt[0] = t*y[0] + y[1]
    dydt[1] = 2*y[0]
    dydt[2] =  2*y[1]+2*t*y[2]
    return dydt

def ode3(t,y) :
    dydt = np.zeros(2)
    dydt[0] = y[1]
    dydt[1] = (5-y[1]**2)*y[1]-y[0]
    return dydt

def ode4(t,y) :
    dydt = -2*y + 2*np.cos(t)*np.sin(2*t)
    return dydt

######################################### my_ODE45

#A=1
#B=2
#arg = A,B
#opts = My.Options()
#sol = My.myODE45(odefcn, np.array([0,5]), np.array([0,0.01]),opts,arg) 

#opts = My.Options()
#opts.odeset('RelTol',1e-10)
#opts.odeset('AbsTol',[1e-10])
#sol = My.myODE45(eq2, np.array([0,5]), np.array([2]),opts)

tspan = [0,10]
y0 = [0]
sol = ode45.ode45(ode1,tspan,y0)

##### real evaluation

#real_sol = np.zeros((1,sol.t.size))
#rel_erreur = np.zeros((1,sol.t.size))
#rel_erreur_norm = np.zeros(sol.t.size)
#for i in range(sol.x.size):
#    real_sol[:,i] = sol2(sol.t[i])
#    rel_erreur[:,i] = np.abs((real_sol[:,i] - sol.y[:,i])/real_sol[:,i])
#    rel_erreur_norm[i] = np.linalg.norm(rel_erreur[:,i])
######################################### plot
    
fig = plt.figure()
plt.title('Basique example')
plt.xlabel('t')
plt.ylabel('y')
plt.plot(sol.t,sol.y[0])
#plt.plot(sol.t,sol.y[1],label='y_2')
plt.show()

######################################### comparaison


#np.savetxt('python_t.txt', sol.t, fmt='%.13e')
#Boolean_t = Compare_data.compare_data('python_t.txt','Matlab_t.txt')
#y = np.transpose(sol.y)
#
#np.savetxt('python_y.txt', y, fmt='%.13e')
#Boolean  = Compare_data.compare_data('python_y.txt','Matlab_y.txt')

#########################################