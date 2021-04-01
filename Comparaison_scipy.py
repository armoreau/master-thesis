import numpy as np
import matplotlib.pyplot as plt
import time
from Options import Options
from scipy.integrate import solve_ivp
from ode45 import ode45

def fun1(t, y):
    return np.array([-0.5 * y])

def sol_fun1(t) : #condition init y[0] = 2
    return np.array([2*np.exp((-1/2)*t)])

def fun2(t,y) :
    return np.cos(t)

def sol_fun2(t) : #condition init y[0] = 0
    return np.sin(t)

def fun3(t, y):
    return np.array([y[1]+y[2], -y[0]+2*y[1]+y[2], y[0]+y[2]])

def sol_fun3(t) : #CI : y[0](0) = 2, y[1](0) = 4, y[2](0) = 6,
    return np.array([-3+5*np.exp(2*t), -3+2*np.exp(t)+5*np.exp(2*t), 3-2*np.exp(t)+5*np.exp(2*t)])

def fun4(t,y) :
    return np.sin(2*t) - np.tan(t)*y

def sol_fun4(t) : #CI : y[0](0) = 1
    return -2*np.cos(t)**2 + 3*np.cos(t)

def fun5(t,y) :
    return (t**2 - t + 1 - t*y)/(t+1)

def sol_fun5(t) : #CI : y[0](1) = 1
    return (t-2) + np.exp(1)*(t+1)*np.exp(-t)

def fun6(t,y) :
    return 4+2*y

def sol_fun6(t) : #CI : y[0](0) = 0
    return 2*np.exp(2*t)-2

###########
length_tspan = 100 # 3 10 100 1000 10000 #base : 10
reltol = 1e-1 #1e-1 1e-3 1e-6 1e-12#base : 1e-3
abstol = 1e-2 #1e-2 1e-6 1e-8 1e-13#base : 1e-6

fig = plt.figure()

################################################################## TEST 1
########### initial condition
t0 = 0
tfinal = 10
tspan = np.linspace(t0,tfinal,length_tspan)
y0 = np.array([2])
myOptions = Options()
myOptions.odeset('RelTol',np.array([reltol]))
myOptions.odeset('AbsTol',abstol*np.ones(len(y0)))
########### scipy aprox
sc_tps1 = time.perf_counter()
sc_sol = solve_ivp(fun1,[t0,tfinal],y0,rtol=reltol,atol=abstol,t_eval=tspan)
sc_tps2 = time.perf_counter()
sc_total_time = sc_tps2-sc_tps1
########### Ode45 aprox
tps1 = time.perf_counter()
sol = ode45(fun1,tspan,y0,myOptions)
tps2 = time.perf_counter()
total_time = tps2-tps1
########### True sol
true_sol = np.zeros([len(y0),len(tspan)])
for i in range(len(tspan)) :
    true_sol[:,i] = sol_fun1(tspan[i])
########### Compute erreur    
abs_error = np.abs(true_sol-sol.y)
#rel_error = abs_error / np.abs(true_sol)
sc_abs_error = np.abs(true_sol-sc_sol.y)
#sc_rel_error = sc_abs_error / np.abs(true_sol)
########### Comparaison Error and time
abs_error_mean = np.mean(abs_error)
sc_abs_error_mean = np.mean(sc_abs_error)
precision_test1 = sc_abs_error_mean/abs_error_mean
speed_test1 = total_time/sc_total_time
########### Plot rel error
#fig = plt.figure()
#plt.title('Test 1 - Comparaison des erreurs absolues')
#plt.xlabel('t')
#plt.ylabel('Erreur absolue')
#plt.plot(sol.t,abs_error[0], label="ode45")
#plt.plot(sc_sol.t,sc_abs_error[0],label="solve_ivp")
#plt.yscale('log')
#plt.legend()

#ax = fig.add_subplot(2, 3, 1)
#ax.plot(sol.t,abs_error[0], label="ode45")
#ax.plot(sc_sol.t,sc_abs_error[0],label="solve_ivp")
#ax.set_yscale('log')
#ax.legend()
################################################################## TEST 2
########### initial condition
t0 = 0
tfinal = 10
tspan = np.linspace(t0,tfinal,length_tspan)
y0 = np.array([0])
myOptions = Options()
myOptions.odeset('RelTol',np.array([reltol]))
myOptions.odeset('AbsTol',abstol*np.ones(len(y0)))
########### scipy aprox
sc_tps1 = time.perf_counter()
sc_sol = solve_ivp(fun2,[t0,tfinal],y0,rtol=reltol,atol=abstol,t_eval=tspan)
sc_tps2 = time.perf_counter()
sc_total_time = sc_tps2-sc_tps1
########### Ode45 aprox
tps1 = time.perf_counter()
sol = ode45(fun2,tspan,y0,myOptions)
tps2 = time.perf_counter()
total_time = tps2-tps1
########### True sol
true_sol = np.zeros([len(y0),len(tspan)])
for i in range(len(tspan)) :
    true_sol[:,i] = sol_fun2(tspan[i])
########### Compute erreur    
abs_error = np.abs(true_sol-sol.y)
#rel_error = abs_error / np.abs(true_sol)
sc_abs_error = np.abs(true_sol-sc_sol.y)
#sc_rel_error = sc_abs_error / np.abs(true_sol)
########### Comparaison Error and time
abs_error_mean = np.mean(abs_error)
sc_abs_error_mean = np.mean(sc_abs_error)
precision_test2 = sc_abs_error_mean/abs_error_mean
speed_test2 = total_time/sc_total_time
########### Plot rel error
#fig = plt.figure()
#plt.title('Test 2 - Comparaison des erreurs absolues')
#plt.xlabel('t')
#plt.ylabel('Erreur absolue')
#plt.plot(sol.t,abs_error[0], label="ode45")
#plt.plot(sc_sol.t,sc_abs_error[0],label="solve_ivp")
#plt.yscale('log')
#plt.legend()

#ax = fig.add_subplot(2, 3, 2)
#ax.plot(sol.t,abs_error[0], label="ode45")
#ax.plot(sc_sol.t,sc_abs_error[0],label="solve_ivp")
#ax.set_yscale('log')
#ax.legend()
################################################################## TEST 3
########### initial condition
t0 = 0
tfinal = 2
tspan = np.linspace(t0,tfinal,length_tspan)
y0 = np.array([2,4,6])
myOptions = Options()
myOptions.odeset('RelTol',np.array([reltol]))
myOptions.odeset('AbsTol',abstol*np.ones(len(y0)))
########### scipy aprox
sc_tps1 = time.perf_counter()
sc_sol = solve_ivp(fun3,[t0,tfinal],y0,rtol=reltol,atol=abstol,t_eval=tspan)
sc_tps2 = time.perf_counter()
sc_total_time = sc_tps2-sc_tps1
########### Ode45 aprox
tps1 = time.perf_counter()
sol = ode45(fun3,tspan,y0,myOptions)
tps2 = time.perf_counter()
total_time = tps2-tps1
########### True sol
true_sol = np.zeros([len(y0),len(tspan)])
for i in range(len(tspan)) :
    true_sol[:,i] = sol_fun3(tspan[i])
########### Compute erreur    
abs_error = np.abs(true_sol-sol.y)
#rel_error = abs_error / np.abs(true_sol)
sc_abs_error = np.abs(true_sol-sc_sol.y)
#sc_rel_error = sc_abs_error / np.abs(true_sol)
########### Comparaison Error and time
abs_error_mean = np.mean(abs_error)
sc_abs_error_mean = np.mean(sc_abs_error)
precision_test3 = sc_abs_error_mean/abs_error_mean
speed_test3 = total_time/sc_total_time
########### Plot rel error
#fig = plt.figure()
#plt.title('Test 3 - Comparaison des erreurs absolues')
#plt.xlabel('t')
#plt.ylabel('Erreur absolue')
#plt.plot(sol.t,abs_error[0], label="ode45")
#plt.plot(sc_sol.t,sc_abs_error[0],label="solve_ivp")
#plt.yscale('log')
#plt.legend()

#ax = fig.add_subplot(2, 3, 3)
#ax.plot(sol.t,abs_error[0], label="ode45")
#ax.plot(sc_sol.t,sc_abs_error[0],label="solve_ivp")
#ax.set_yscale('log')
#ax.legend()


ax = fig.add_subplot(1, 2, 1)
ax.plot(sol.t,abs_error[0], label="ode45")
ax.plot(sc_sol.t,sc_abs_error[0],label="solve_ivp")
ax.set_yscale('log')
ax.set_title('Précision Basse')
ax.legend()

reltol = 1e-3 #1e-1 1e-3 1e-6 1e-12#base : 1e-3
abstol = 1e-6 #1e-2 1e-6 1e-8 1e-13#base : 1e-6

myOptions = Options()
myOptions.odeset('RelTol',np.array([reltol]))
myOptions.odeset('AbsTol',abstol*np.ones(len(y0)))
########### scipy aprox
sc_tps1 = time.perf_counter()
sc_sol = solve_ivp(fun3,[t0,tfinal],y0,rtol=reltol,atol=abstol,t_eval=tspan)
sc_tps2 = time.perf_counter()
sc_total_time = sc_tps2-sc_tps1
########### Ode45 aprox
tps1 = time.perf_counter()
sol = ode45(fun3,tspan,y0,myOptions)
tps2 = time.perf_counter()
total_time = tps2-tps1
########### True sol
true_sol = np.zeros([len(y0),len(tspan)])
for i in range(len(tspan)) :
    true_sol[:,i] = sol_fun3(tspan[i])
########### Compute erreur    
abs_error = np.abs(true_sol-sol.y)
#rel_error = abs_error / np.abs(true_sol)
sc_abs_error = np.abs(true_sol-sc_sol.y)
#sc_rel_error = sc_abs_error / np.abs(true_sol)
########### Comparaison Error and time
abs_error_mean = np.mean(abs_error)
sc_abs_error_mean = np.mean(sc_abs_error)
precision_test3 = sc_abs_error_mean/abs_error_mean
speed_test3 = total_time/sc_total_time


ax = fig.add_subplot(1, 2, 2)
ax.plot(sol.t,abs_error[0], label="ode45")
ax.plot(sc_sol.t,sc_abs_error[0],label="solve_ivp")
ax.set_yscale('log')
ax.set_title('Précision Moyenne')
ax.legend()

reltol = 1e-6 #1e-1 1e-3 1e-6 1e-12#base : 1e-3
abstol = 1e-8 #1e-2 1e-6 1e-8 1e-13#base : 1e-6

myOptions = Options()
myOptions.odeset('RelTol',np.array([reltol]))
myOptions.odeset('AbsTol',abstol*np.ones(len(y0)))
########### scipy aprox
sc_tps1 = time.perf_counter()
sc_sol = solve_ivp(fun3,[t0,tfinal],y0,rtol=reltol,atol=abstol,t_eval=tspan)
sc_tps2 = time.perf_counter()
sc_total_time = sc_tps2-sc_tps1
########### Ode45 aprox
tps1 = time.perf_counter()
sol = ode45(fun3,tspan,y0,myOptions)
tps2 = time.perf_counter()
total_time = tps2-tps1
########### True sol
true_sol = np.zeros([len(y0),len(tspan)])
for i in range(len(tspan)) :
    true_sol[:,i] = sol_fun3(tspan[i])
########### Compute erreur    
abs_error = np.abs(true_sol-sol.y)
#rel_error = abs_error / np.abs(true_sol)
sc_abs_error = np.abs(true_sol-sc_sol.y)
#sc_rel_error = sc_abs_error / np.abs(true_sol)
########### Comparaison Error and time
abs_error_mean = np.mean(abs_error)
sc_abs_error_mean = np.mean(sc_abs_error)
precision_test3 = sc_abs_error_mean/abs_error_mean
speed_test3 = total_time/sc_total_time

fig = plt.figure()

ax = fig.add_subplot(1, 2, 1)
ax.plot(sol.t,abs_error[0], label="ode45")
ax.plot(sc_sol.t,sc_abs_error[0],label="solve_ivp")
ax.set_yscale('log')
ax.set_title('Précision Haute')
ax.legend()

reltol = 1e-12 #1e-1 1e-3 1e-6 1e-12#base : 1e-3
abstol = 1e-13 #1e-2 1e-6 1e-8 1e-13#base : 1e-6

myOptions = Options()
myOptions.odeset('RelTol',np.array([reltol]))
myOptions.odeset('AbsTol',abstol*np.ones(len(y0)))
########### scipy aprox
sc_tps1 = time.perf_counter()
sc_sol = solve_ivp(fun3,[t0,tfinal],y0,rtol=reltol,atol=abstol,t_eval=tspan)
sc_tps2 = time.perf_counter()
sc_total_time = sc_tps2-sc_tps1
########### Ode45 aprox
tps1 = time.perf_counter()
sol = ode45(fun3,tspan,y0,myOptions)
tps2 = time.perf_counter()
total_time = tps2-tps1
########### True sol
true_sol = np.zeros([len(y0),len(tspan)])
for i in range(len(tspan)) :
    true_sol[:,i] = sol_fun3(tspan[i])
########### Compute erreur    
abs_error = np.abs(true_sol-sol.y)
#rel_error = abs_error / np.abs(true_sol)
sc_abs_error = np.abs(true_sol-sc_sol.y)
#sc_rel_error = sc_abs_error / np.abs(true_sol)
########### Comparaison Error and time
abs_error_mean = np.mean(abs_error)
sc_abs_error_mean = np.mean(sc_abs_error)
precision_test3 = sc_abs_error_mean/abs_error_mean
speed_test3 = total_time/sc_total_time

ax = fig.add_subplot(1, 2, 2)
ax.plot(sol.t,abs_error[0], label="ode45")
ax.plot(sc_sol.t,sc_abs_error[0],label="solve_ivp")
ax.set_yscale('log')
ax.set_title('Précision très Haute')
ax.legend()

################################################################## TEST 4
########### initial condition
t0 = 0
tfinal = np.pi/2 - 0.3
tspan = np.linspace(t0,tfinal,length_tspan)
y0 = np.array([1])
myOptions = Options()
myOptions.odeset('RelTol',np.array([reltol]))
myOptions.odeset('AbsTol',abstol*np.ones(len(y0)))
########### scipy aprox
sc_tps1 = time.perf_counter()
sc_sol = solve_ivp(fun4,[t0,tfinal],y0,rtol=reltol,atol=abstol,t_eval=tspan)
sc_tps2 = time.perf_counter()
sc_total_time = sc_tps2-sc_tps1
########### Ode45 aprox
tps1 = time.perf_counter()
sol = ode45(fun4,tspan,y0,myOptions)
tps2 = time.perf_counter()
total_time = tps2-tps1
########### True sol
true_sol = np.zeros([len(y0),len(tspan)])
for i in range(len(tspan)) :
    true_sol[:,i] = sol_fun4(tspan[i])
########### Compute erreur    
abs_error = np.abs(true_sol-sol.y)
#rel_error = abs_error / np.abs(true_sol)
sc_abs_error = np.abs(true_sol-sc_sol.y)
#sc_rel_error = sc_abs_error / np.abs(true_sol)
########### Comparaison Error and time
abs_error_mean = np.mean(abs_error)
sc_abs_error_mean = np.mean(sc_abs_error)
precision_test4 = sc_abs_error_mean/abs_error_mean
speed_test4 = total_time/sc_total_time
########### Plot rel error
#fig = plt.figure()
#plt.title('Test 4 - Comparaison des erreurs absolues')
#plt.xlabel('t')
#plt.ylabel('Erreur absolue')
#plt.plot(sol.t,abs_error[0], label="ode45")
#plt.plot(sc_sol.t,sc_abs_error[0],label="solve_ivp")
#plt.yscale('log')
#plt.legend()

#ax = fig.add_subplot(2, 3, 4)
#ax.plot(sol.t,abs_error[0], label="ode45")
#ax.plot(sc_sol.t,sc_abs_error[0],label="solve_ivp")
#ax.set_yscale('log')
#ax.legend()
################################################################## TEST 5
########### initial condition
t0 = 1
tfinal = 10
tspan = np.linspace(t0,tfinal,length_tspan)
y0 = np.array([1])
myOptions = Options()
myOptions.odeset('RelTol',np.array([reltol]))
myOptions.odeset('AbsTol',abstol*np.ones(len(y0)))
########### scipy aprox
sc_tps1 = time.perf_counter()
sc_sol = solve_ivp(fun5,[t0,tfinal],y0,rtol=reltol,atol=abstol,t_eval=tspan)
sc_tps2 = time.perf_counter()
sc_total_time = sc_tps2-sc_tps1
########### Ode45 aprox
tps1 = time.perf_counter()
sol = ode45(fun5,tspan,y0,myOptions)
tps2 = time.perf_counter()
total_time = tps2-tps1
########### True sol
true_sol = np.zeros([len(y0),len(tspan)])
for i in range(len(tspan)) :
    true_sol[:,i] = sol_fun5(tspan[i])
########### Compute erreur    
abs_error = np.abs(true_sol-sol.y)
#rel_error = abs_error / np.abs(true_sol)
sc_abs_error = np.abs(true_sol-sc_sol.y)
#sc_rel_error = sc_abs_error / np.abs(true_sol)
########### Comparaison Error and time
abs_error_mean = np.mean(abs_error)
sc_abs_error_mean = np.mean(sc_abs_error)
precision_test5 = sc_abs_error_mean/abs_error_mean
speed_test5 = total_time/sc_total_time
########### Plot rel error
#fig = plt.figure()
#plt.title('Test 5 - Comparaison des erreurs absolues')
#plt.xlabel('t')
#plt.ylabel('Erreur absolue')
#plt.plot(sol.t,abs_error[0], label="ode45")
#plt.plot(sc_sol.t,sc_abs_error[0],label="solve_ivp")
#plt.yscale('log')
#plt.legend()

#ax = fig.add_subplot(2, 3, 5)
#ax.plot(sol.t,abs_error[0], label="ode45")
#ax.plot(sc_sol.t,sc_abs_error[0],label="solve_ivp")
#ax.set_yscale('log')
#ax.legend()
################################################################## TEST 6
########### initial condition
t0 = 0
tfinal = 2
tspan = np.linspace(t0,tfinal,length_tspan)
y0 = np.array([0])
myOptions = Options()
myOptions.odeset('RelTol',np.array([reltol]))
myOptions.odeset('AbsTol',abstol*np.ones(len(y0)))
########### scipy aprox
sc_tps1 = time.perf_counter()
sc_sol = solve_ivp(fun6,[t0,tfinal],y0,rtol=reltol,atol=abstol,t_eval=tspan)
sc_tps2 = time.perf_counter()
sc_total_time = sc_tps2-sc_tps1
########### Ode45 aprox
tps1 = time.perf_counter()
sol = ode45(fun6,tspan,y0,myOptions)
tps2 = time.perf_counter()
total_time = tps2-tps1
########### True sol
true_sol = np.zeros([len(y0),len(tspan)])
for i in range(len(tspan)) :
    true_sol[:,i] = sol_fun6(tspan[i])
########### Compute erreur    
abs_error = np.abs(true_sol-sol.y)
#rel_error = abs_error / np.abs(true_sol)
sc_abs_error = np.abs(true_sol-sc_sol.y)
#sc_rel_error = sc_abs_error / np.abs(true_sol)
########### Comparaison Error and time
abs_error_mean = np.mean(abs_error)
sc_abs_error_mean = np.mean(sc_abs_error)
precision_test6 = sc_abs_error_mean/abs_error_mean
speed_test6 = total_time/sc_total_time
########### Plot rel error
#fig = plt.figure()
#plt.title('Test 6 - Comparaison des erreurs absolues')
#plt.xlabel('t')
#plt.ylabel('Erreur absolue')
#plt.plot(sol.t,abs_error[0], label="ode45")
#plt.plot(sc_sol.t,sc_abs_error[0],label="solve_ivp")
#plt.yscale('log')
#plt.legend()

#ax = fig.add_subplot(2, 3, 6)
#ax.plot(sol.t,abs_error[0], label="ode45")
#ax.plot(sc_sol.t,sc_abs_error[0],label="solve_ivp")
#ax.set_yscale('log')
#ax.legend()

plt.show()

##########################
speed = np.array([speed_test1, speed_test2, speed_test3, speed_test4, speed_test5, speed_test6])
precision = np.array([precision_test1, precision_test2, precision_test3, precision_test4, precision_test5, precision_test6])

