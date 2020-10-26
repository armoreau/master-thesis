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

######################################### solve_ivp 

tps1 = time.perf_counter()
solu = solve_ivp(eq2, [0, 100], [2])
tps2 = time.perf_counter()
total_time = tps2 - tps1

######################################### my_ODE45

my_tps1 = time.perf_counter()
[t, y] = myODE45(eq2, [0, 100], [2], None)
my_tps2 = time.perf_counter()
my_total_time = my_tps2 - my_tps1

######################################### real evaluation

real_sol = np.zeros((np.size(t)))

for i in range(np.size(t)) :
    real_sol[i] = sol2(t[i])
    
######################################### error

error = np.abs(real_sol - solu.y)
my_error = np.abs(real_sol - y)

#########################################

print('total_time_solve_ivp = ')
print(total_time)
print('my_total_time_myODE45 = ')
print(my_total_time)

#########################################
fig = plt.figure()
plt.title('EDO : dy/dt = -0.5 * y         y(0) = 2         sol : y(t) = 2*exp(-t/2)')
plt.xlabel('t')
plt.ylabel('approximation_error')
plt.yscale('log')
plt.plot(t,error[0],label='solve_ivp_error')
plt.plot(t,my_error[0],label='myODE45_error')
plt.legend()
plt.show()