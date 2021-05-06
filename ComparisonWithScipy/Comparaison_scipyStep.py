import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.integrate import solve_ivp

#Add parent folder to the path. Code taken from https://codeolives.com/2020/01/10/python-reference-module-in-parent-directory/
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from odeoptions import Odeoptions
from ode45_scipyStep import ode45_scipyStep

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
reltol = 1e-3 #1e-1 1e-3 1e-6 1e-12#base : 1e-3
abstol = 1e-6 #1e-2 1e-6 1e-8 1e-13#base : 1e-6
myOptions = Odeoptions()
myOptions.odeset('Refine',1)
myOptions.odeset('RelTol',np.array([reltol]))
################################################################## TEST 1
########### initial condition
t0 = 0
tfinal = 100
tspan = [t0,tfinal]
y0 = np.array([2])
myOptions.odeset('AbsTol',abstol*np.ones(len(y0)))
########### scipy aprox
sc_tps1 = time.perf_counter()
sc_sol = solve_ivp(fun1,[t0,tfinal],y0,rtol=reltol,atol=abstol)
sc_tps2 = time.perf_counter()
sc_total_time = sc_tps2-sc_tps1

scipyStep = np.zeros(len(sc_sol.t)-1)
for i in range(len(sc_sol.t)-1) :
    scipyStep[i] = sc_sol.t[i+1] - sc_sol.t[i]
########### Ode45 aprox
tps1 = time.perf_counter()
sol = ode45_scipyStep(fun1,tspan,y0,scipyStep,myOptions)
tps2 = time.perf_counter()
total_time = tps2-tps1
###########
speed_test1 = sc_total_time/total_time
##################################### Assure Equal precision
########### True sol
true_sol = np.zeros([len(y0),len(sol.t)])
for i in range(len(sol.t)) :
    true_sol[:,i] = sol_fun1(sol.t[i])
########### Compute erreur    
abs_error = np.abs(true_sol-sol.y)
sc_abs_error = np.abs(true_sol-sc_sol.y)
########### Comparaison Error and time
abs_error_mean = np.mean(abs_error)
sc_abs_error_mean = np.mean(sc_abs_error)
precision_test1 = sc_abs_error_mean/abs_error_mean
################################################################## TEST 2
########### initial condition
t0 = 0
tfinal = 10
tspan = [t0,tfinal]
y0 = np.array([0])
myOptions.odeset('AbsTol',abstol*np.ones(len(y0)))
########### scipy aprox
sc_tps1 = time.perf_counter()
sc_sol = solve_ivp(fun2,[t0,tfinal],y0,rtol=reltol,atol=abstol)
sc_tps2 = time.perf_counter()
sc_total_time = sc_tps2-sc_tps1

scipyStep = np.zeros(len(sc_sol.t)-1)
for i in range(len(sc_sol.t)-1) :
    scipyStep[i] = sc_sol.t[i+1] - sc_sol.t[i]
########### Ode45 aprox
tps1 = time.perf_counter()
sol = ode45_scipyStep(fun2,tspan,y0,scipyStep,myOptions)
tps2 = time.perf_counter()
total_time = tps2-tps1
###########
speed_test2 = sc_total_time/total_time
################################################################## TEST 3
########### initial condition
t0 = 0
tfinal = 2
tspan = [t0,tfinal]
y0 = np.array([2,4,6])
myOptions.odeset('AbsTol',abstol*np.ones(len(y0)))
########### scipy aprox
sc_tps1 = time.perf_counter()
sc_sol = solve_ivp(fun3,[t0,tfinal],y0,rtol=reltol,atol=abstol)
sc_tps2 = time.perf_counter()
sc_total_time = sc_tps2-sc_tps1

scipyStep = np.zeros(len(sc_sol.t)-1)
for i in range(len(sc_sol.t)-1) :
    scipyStep[i] = sc_sol.t[i+1] - sc_sol.t[i]
########### Ode45 aprox
tps1 = time.perf_counter()
sol = ode45_scipyStep(fun3,tspan,y0,scipyStep,myOptions)
tps2 = time.perf_counter()
total_time = tps2-tps1
###########
speed_test3 = sc_total_time/total_time
################################################################## TEST 4
########### initial condition
t0 = 0
tfinal = np.pi/2 - 0.3
tspan = [t0,tfinal]
y0 = np.array([1])
myOptions.odeset('AbsTol',abstol*np.ones(len(y0)))
########### scipy aprox
sc_tps1 = time.perf_counter()
sc_sol = solve_ivp(fun2,[t0,tfinal],y0,rtol=reltol,atol=abstol)
sc_tps2 = time.perf_counter()
sc_total_time = sc_tps2-sc_tps1

scipyStep = np.zeros(len(sc_sol.t)-1)
for i in range(len(sc_sol.t)-1) :
    scipyStep[i] = sc_sol.t[i+1] - sc_sol.t[i]
########### Ode45 aprox
tps1 = time.perf_counter()
sol = ode45_scipyStep(fun2,tspan,y0,scipyStep,myOptions)
tps2 = time.perf_counter()
total_time = tps2-tps1
###########
speed_test4 = sc_total_time/total_time
################################################################## TEST 5
########### initial condition
t0 = 1
tfinal = 10
tspan = [t0,tfinal]
y0 = np.array([1])
myOptions.odeset('AbsTol',abstol*np.ones(len(y0)))
########### scipy aprox
sc_tps1 = time.perf_counter()
sc_sol = solve_ivp(fun2,[t0,tfinal],y0,rtol=reltol,atol=abstol)
sc_tps2 = time.perf_counter()
sc_total_time = sc_tps2-sc_tps1

scipyStep = np.zeros(len(sc_sol.t)-1)
for i in range(len(sc_sol.t)-1) :
    scipyStep[i] = sc_sol.t[i+1] - sc_sol.t[i]
########### Ode45 aprox
tps1 = time.perf_counter()
sol = ode45_scipyStep(fun2,tspan,y0,scipyStep,myOptions)
tps2 = time.perf_counter()
total_time = tps2-tps1
###########
speed_test5 = sc_total_time/total_time
################################################################## TEST 6
########### initial condition
t0 = 0
tfinal = 2
tspan = [t0,tfinal]
y0 = np.array([0])
myOptions.odeset('AbsTol',abstol*np.ones(len(y0)))
########### scipy aprox
sc_tps1 = time.perf_counter()
sc_sol = solve_ivp(fun2,[t0,tfinal],y0,rtol=reltol,atol=abstol)
sc_tps2 = time.perf_counter()
sc_total_time = sc_tps2-sc_tps1

scipyStep = np.zeros(len(sc_sol.t)-1)
for i in range(len(sc_sol.t)-1) :
    scipyStep[i] = sc_sol.t[i+1] - sc_sol.t[i]
########### Ode45 aprox
tps1 = time.perf_counter()
sol = ode45_scipyStep(fun2,tspan,y0,scipyStep,myOptions)
tps2 = time.perf_counter()
total_time = tps2-tps1
###########
speed_test6 = sc_total_time/total_time


################################################
speed = np.array([speed_test1, speed_test2, speed_test3, speed_test4, speed_test5, speed_test6])