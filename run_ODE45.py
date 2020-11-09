from scipy.integrate import solve_ivp
from MyODE45 import myODE45
import numpy as np
import matplotlib.pyplot as plt
import time


def eq1(t, y):
    return [y[1]+y[2], -y[0]+2*y[1]+y[2], y[0]+y[2]]
def sol1(t) : #CI : y[0](0) = 2, y[1](0) = 4, y[2](0) = 6,
    return [-3+5*np.exp(2*t), -3+2*np.exp(t)+5*np.exp(2*t), 3-2*np.exp(t)+5*np.exp(2*t)]

def eq2(t, y):
    return -0.5 * y[0]
def sol2(t) : #condition init y[0] = 2
    return 2*np.exp((-1/2)*t) 


######################################### my_ODE45

#t_eval = np.linspace(0,10,1000)
my_tps1 = time.perf_counter()
[t, y] = myODE45(eq1, np.array([0, 100]), np.array([2,4,6])) #np.array([0, 100])
my_tps2 = time.perf_counter()
my_total_time = my_tps2 - my_tps1

##### real evaluation

real_sol = np.zeros((3,t.size))
rel_erreur = np.zeros((3,t.size))
rel_erreur_norm = np.zeros(t.size)
for i in range(t.size):
    real_sol[:,i] = sol1(t[i])
    rel_erreur[:,i] = np.abs((real_sol[:,i] - y[:,i])/real_sol[:,i])
    rel_erreur_norm[i] = np.linalg.norm(rel_erreur[:,i])


######################################### print

#print(rel_erreur_norm)
print(t)
print(y)
print(my_total_time)

#########################################
fig = plt.figure()
plt.title('Relative error')
plt.xlabel('t')
plt.ylabel('relative error')
plt.yscale('log')
plt.plot(t,rel_erreur_norm,label='x')
plt.legend()

fig = plt.figure()
plt.xlabel('t')
plt.ylabel('Y_out(t)')
plt.plot(t,y[0,:],label='x(t)')
plt.plot(t,y[1,:],label='y(t)')
plt.plot(t,y[2,:],label='z(t)')
plt.legend()

plt.show()