from scipy.integrate import solve_ivp
from MyODE45 import myODE45
import numpy as np
import matplotlib.pyplot as plt
import time


def diff_expo(t, y): return -0.5 * y
def expo(t) : return 2*np.exp((-1/2)*t)

def diff_sin(t,y) : return np.cos(t)
def sin(t) : return np.sin(t)

Nbr_pts = 101
t = np.linspace(0, 10, Nbr_pts)

######################################### Evaluation exacte
true_sol_expo = np.zeros(Nbr_pts)
true_sol_sin = np.zeros(Nbr_pts)
for i in range(Nbr_pts) :
    true_sol_expo[i] = expo(t[i])
    true_sol_sin[i] = sin(t[i])
    
    
######################################### Approximation par solve_ivp
print("===================")

tps1_expo = time.perf_counter()
sol_expo = solve_ivp(diff_expo, [0, 10], [2], t_eval=t)
tps2_expo = time.perf_counter()

tps1_sin = time.perf_counter()
sol_sin = solve_ivp(diff_sin, [0, 10], [0], t_eval=t)
tps2_sin = time.perf_counter()

total_time_expo = tps2_expo - tps1_expo
total_time_sin = tps2_sin - tps1_sin 

erreur_expo = np.abs(true_sol_expo - sol_expo.y[0])
erreur_sin = np.abs(true_sol_sin - sol_sin.y[0])
    
######################################### Mon approximation
print("===================")

my_tps1_expo = time.perf_counter()
(my_sol_expo_t,my_sol_expo_y) = myODE45(diff_expo, [0, 10], [2], t)
my_tps2_expo = time.perf_counter()

my_tps1_sin = time.perf_counter()
(my_sol_sin_t,my_sol_sin_y) = myODE45(diff_sin, [0, 10], [0], t)
my_tps2_sin = time.perf_counter()

my_total_time_expo = my_tps2_expo - my_tps1_expo
my_total_time_sin = my_tps2_sin - my_tps1_sin 

my_erreur_expo = np.abs(true_sol_expo - my_sol_expo_y[0])
my_erreur_sin = np.abs(true_sol_sin - my_sol_sin_y[0])


######################################### Plot

fig = plt.figure()
plt.title('EDO : dy/dt = -0.5 * y         y(0) = 2         sol : y(t) = 2*exp(-t/2)')
plt.xlabel('t')
plt.ylabel('approximation_error')
plt.yscale('log')
plt.plot(t,erreur_expo,label='solve_ivp_error')
plt.plot(t,my_erreur_expo,label='myODE45_error')
plt.legend()

fig2 = plt.figure()
plt.title('EDO : dy/dt = cos(t)       y(0) = 0        sol : y(t) = sin(t)')
plt.xlabel('t')
plt.ylabel('approximation_error')
plt.yscale('log')
plt.plot(t,erreur_sin,label='solve_ivp_error')
plt.plot(t,my_erreur_sin,label='myODE45_error')
plt.legend()

#plt.show()


print('total_time_expo = ')
print(total_time_expo)
print('my_total_time_expo = ')
print(my_total_time_expo)

print('total_time_sin = ')
print(total_time_sin)
print('my_total_time_sin = ')
print(my_total_time_sin)

######################################### Nbrs laisser choisit par le solveur

######################################### Approximation par solve_ivp
print("===================")

tps1_expo = time.perf_counter()
sol_expo = solve_ivp(diff_expo, [0, 10], [2])
tps2_expo = time.perf_counter()

tps1_sin = time.perf_counter()
sol_sin = solve_ivp(diff_sin, [0, 10], [0])
tps2_sin = time.perf_counter()

total_time_expo = tps2_expo - tps1_expo
total_time_sin = tps2_sin - tps1_sin

############ Evaluation exacte
true_sol_expo = np.zeros(np.size(sol_expo.t))
true_sol_sin = np.zeros(np.size(sol_sin.t))

for i in range(np.size(sol_expo.t)) :
    true_sol_expo[i] = expo(sol_expo.t[i])
for i in range(np.size(sol_sin.t)) :
    true_sol_sin[i] = sin(sol_sin.t[i])
############

erreur_expo = np.abs(true_sol_expo - sol_expo.y[0])
erreur_sin = np.abs(true_sol_sin - sol_sin.y[0])
    
######################################### Mon approximation
print("===================")

my_tps1_expo = time.perf_counter()
(my_sol_expo_t,my_sol_expo_y) = myODE45(diff_expo, [0, 10], [2], None)
my_tps2_expo = time.perf_counter()

my_tps1_sin = time.perf_counter()
(my_sol_sin_t,my_sol_sin_y) = myODE45(diff_sin, [0, 10], [0], None)
my_tps2_sin = time.perf_counter()

my_total_time_expo = my_tps2_expo - my_tps1_expo
my_total_time_sin = my_tps2_sin - my_tps1_sin

############ Evaluation exacte
true_sol_expo = np.zeros(np.size(my_sol_expo_t))
true_sol_sin = np.zeros(np.size(my_sol_sin_t))
for i in range(np.size(my_sol_expo_t)) :
    true_sol_expo[i] = expo(my_sol_expo_t[0][i])
for i in range(np.size(my_sol_sin_t)) :
    true_sol_sin[i] = sin(my_sol_sin_t[0][i])
############

my_erreur_expo = np.abs(true_sol_expo - my_sol_expo_y[0])
my_erreur_sin = np.abs(true_sol_sin - my_sol_sin_y[0])

fig3 = plt.figure()
plt.title('EDO : dy/dt = -0.5 * y         y(0) = 2         sol : y(t) = 2*exp(-t/2)')
plt.xlabel('t')
plt.ylabel('approximation_error')
plt.yscale('log')
plt.plot(sol_expo.t,erreur_expo,label='solve_ivp_error')
plt.plot(my_sol_expo_t[0],my_erreur_expo,label='myODE45_error')
plt.legend()

fig4 = plt.figure()
plt.title('EDO : dy/dt = cos(t)       y(0) = 0        sol : y(t) = sin(t)')
plt.xlabel('t')
plt.ylabel('approximation_error')
plt.yscale('log')
plt.plot(sol_sin.t,erreur_sin,label='solve_ivp_error')
plt.plot(my_sol_sin_t[0],my_erreur_sin,label='myODE45_error')
plt.legend()

print('total_time_expo = ')
print(total_time_expo)
print('my_total_time_expo = ')
print(my_total_time_expo)

print('total_time_sin = ')
print(total_time_sin)
print('my_total_time_sin = ')
print(my_total_time_sin)

plt.show()