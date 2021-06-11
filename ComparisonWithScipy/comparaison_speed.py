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

######################################### function test
### First order fction
### trigonometrique
def f1(t,y) :
    dydt = np.cos(t)
    return dydt

def f2(t,y) :
    dydt = np.sin(t)
    return dydt

def f3(t,y) :
    dydt = np.sin(t) - np.cos(t)
    return dydt

### exponential
def f4(t,y) :
    dydt = np.exp(1)**(t)
    return dydt

def f5(t,y) :
    dydt = y
    return dydt

def f6(t,y) :
    dydt = [-0.5*y]
    return dydt

def f7(t,y) :
    dydt = 0.5*y
    return dydt

### polynomial
def f8(t,y) :
    dydt = 2*t
    return dydt

def f9(t,y) :
    dydt = t**2 + 2*t
    return dydt

def f10(t,y) :
    dydt = 5*t**5 - t**2 + 2*t
    return dydt

### other
def f11(t,y) :
    dydt = y + t
    return dydt

def f12(t,y) :
    dydt = np.sin(t) + t
    return dydt

def f13(t,y) :
    dydt = np.sin(t) + y
    return dydt

def f14(t,y) :
    dydt = 2*y+t**2
    return dydt

def f15(t,y) :
    dydt = np.sin(t)**2
    return dydt

#### second order or more fction

def f16(t,y) :
    dydt = np.zeros(2)
    dydt[0] = y[1]
    dydt[1] = (1-y[0]**2)*y[1]-y[0]
    return dydt

def f17(t,y) :
    dydt = np.zeros(2)
    dydt[0] = t
    dydt[1] = y[1]
    return dydt

def f18(t,y) :
    dydt = np.zeros(2)
    dydt[0] = y[1]
    dydt[1] = y[0]
    return dydt

def f19(t,y) :
    dydt = np.zeros(2)
    dydt[0] = y[1]
    dydt[1] = t*y[0]
    return dydt

def f20(t,y) :
    dydt = np.zeros(2)
    dydt[0] = 2*y[0]
    dydt[1] = 1/(y[0]+100)
    return dydt

def f21(t,y) :
    dydt = np.zeros(3)
    dydt[0] = t
    dydt[1] = y[0]
    dydt[2] = 2*y[1]
    return dydt

def f22(t,y) :
    dydt = np.zeros(3)
    dydt[0] = t
    dydt[1] = y[0]
    dydt[2] = y[0] + 1
    return dydt

#### additionnal argument fction
def f23(t,y,A,B) :
    dydt = np.zeros(2)
    dydt[0] = y[1]
    dydt[1] = (A/B)*t* y[0]
    return dydt

def f24(t,y,A,B,C) :
    dydt = (A+B)*y/C
    return dydt

def f25(t,y,A,B,C) :
    dydt = np.sin((A+B)*t/C)
    return dydt

#### event fction

def f26(t,y):
    dydt = np.zeros(2)
    dydt[0] = y[1]
    dydt[1] = -9.8
    return dydt

def events_f26(t,y):
    value = [y[0]]
    isterminal = [1]
    direction = [1]
    return [value,isterminal,direction]

def events_f26_scipy(t,y):
    return y[0]

def f27(t,y) :
    mu = 1 / 82.45
    mustar = 1 - mu
    r13 = ((y[0] + mu)**2 + y[1]**2) ** 1.5
    r23 = ((y[0] - mustar)**2 + y[1]**2) ** 1.5
    
    dydt = np.zeros(4)
    dydt[0] = y[2]
    dydt[1] = y[3]
    dydt[2] = 2*y[3] + y[0] - mustar*((y[0]+mu)/r13) - mu*((y[0]-mustar)/r23)
    dydt[3] = -2*y[2] + y[1] - mustar*(y[1]/r13) - mu*(y[1]/r23)
    return dydt

def events_f27(t,y):
    y0 = np.array([1.2, 0, 0, -1.04935750983031990726])
    dDSQdt = 2 * np.dot((y0[0:2]-y0[0:2]),y0[2:4])
    value = [dDSQdt, dDSQdt]
    isterminal = [1, 0]
    direction = [1, 1]
    return [value,isterminal,direction]

def events_f27_scipy(t,y):
    y0 = np.array([1.2, 0, 0, -1.04935750983031990726])
    dDSQdt = 2 * np.dot((y0[0:2]-y0[0:2]),y0[2:4])
    return dDSQdt

def compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C) :
    
    scipy_time = np.zeros([27])
    ode45_time = np.zeros([27])
    i=0
    
    arg = None
    y0 = y0_1
    opts.Refine = 1
    
    #### Test 1
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f1,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
    
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f1, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 2
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f2,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
        
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f2, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 3
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f3,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
        
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f3, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 4
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f4,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
        
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f4, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 5
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f5,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
       
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f5, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 6
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f6,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
         
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f6, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 7
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f7,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
         
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f7, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 8
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f8,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
         
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f8, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 9
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f9,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
         
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f9, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 10
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f10,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
         
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f10, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 11
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f11,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
         
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f11, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 12
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f12,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
         
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f12, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 13
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f13,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
         
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f13, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 14
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f14,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
         
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f14, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 15
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f15,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]    
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f15, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    ### len(y0)=2
    y0 = y0_2

    #### Test 16
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f16,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
         
    t1 = time.perf_counter() 
    sol = ode45_scipyStep(f16, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 17
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f17,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
         
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f17, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 18
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f18,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
         
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f18, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 19
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f19,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
       
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f19, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 20
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f20,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
    
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f20, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    ### len(y0)=3
    y0 = y0_3

    #### Test 21
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f21,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
        
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f21, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 22
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f22,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
        
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f22, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    ### arg fction
    arg = A,B
    y0 = y0_2

    #### Test 23
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f23,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
        
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f23, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1
    
    ### arg fction
    arg = A,B,C
    y0 = y0_1

    #### Test 24
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f24,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
        
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f24, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1

    #### Test 25
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f25,tspan,y0,events=opts.Events,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
        
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f25, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1
    
    arg = None
    
    ### event fction
    y0 = [0.0, 20.0]
    opts.odeset('Events',events_f26)

    #### Test 26
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f26,tspan,y0,events=events_f26_scipy,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
    
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f26, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1
    
    ### event fction
    y0 = [1.2, 0, 0, -1.04935750983031990726]
    opts.odeset('Events',events_f27)

    #### Test 27
    scipy_t1 = time.perf_counter()
    scipy_sol = solve_ivp(f27,tspan,y0,events=events_f27_scipy,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    scipy_t2 = time.perf_counter()
    scipy_time[i] = scipy_t2 - scipy_t1
    scipyStep = np.zeros(len(scipy_sol.t)-1)
    for j in range(len(scipy_sol.t)-1) :
        scipyStep[j] = scipy_sol.t[j+1] - scipy_sol.t[j]
        
    t1 = time.perf_counter()
    sol = ode45_scipyStep(f27, tspan, y0, scipyStep, opts, arg)
    t2 = time.perf_counter()
    ode45_time[i] = t2-t1
    i = i+1
    
    return ode45_time/scipy_time

######################################### perform test with random input
nbr_test = 27
nbr_Input = 35
temps = np.zeros([nbr_Input, nbr_test])

#### INPUT 1
INPUT = 0
tspan = [ 4.60789919, 14.85334451]
y0_1 = [4.53897499]
y0_2 = [1.20051259, 1.17666892]
y0_3 = [ 0.45052893,  7.03086018, -1.50526611]
A = -2.83009697
B = -4.97673659
C = -3.50342477
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 2
INPUT = INPUT + 1
tspan = [ 9.91082815, 11.84808434]
y0_1 = [-6.38784938]
y0_2 = [-1.60628946, -1.18921585]
y0_3 = [-7.09018137,  4.81021195, -7.82701511]
A = 2.25480317
B = -0.94979655
C = 0.1968337
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 3
INPUT = INPUT + 1
tspan =  [-1.41580028, 7.95851231]
y0_1 = [-1.94429498]
y0_2 = [4.96831322, -7.15530037]
y0_3 = [-2.33703264, -2.00682155, -7.14666703]
A = 4.75708095
B = 3.95418388
C = -4.72783249
opts = Odeoptions()
opts.odeset('RelTol',1.3818760883438436e-06)
opts.odeset('AbsTol',2.0635002062377786e-08)
opts.odeset('InitialStep',0.23252681078089357)
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 4
INPUT = INPUT + 1
tspan = [-2.08547886,  2.33452522]
y0_1 = [-3.76904574]
y0_2 = [-5.24019054,  4.65499285]
y0_3 = [-6.51846199,  3.11856858,  1.27300755]
A = 1.30937634
B = -1.17558159
C = -1.51165755
opts = Odeoptions()
opts.odeset('AbsTol',1.4482218440447706e-09)
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 5
INPUT = INPUT + 1
tspan = [-0.32081691, 2.54055559]
y0_1 = [-0.01349065]
y0_2 = [-0.76687766, 3.93978074]
y0_3 = [2.98288271, 3.48166087, 7.2604428 ]
A = 6.61349864e-01
B = -3.65342548e-03
C = -4.36285600e+00
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 6
INPUT = INPUT + 1
tspan = [10.20618698, 14.25437019]
y0_1 = [-2.51324011]
y0_2 = [-6.18446334,  7.34621481]
y0_3 = [-4.52846858,  2.46373529,  2.13908458]
A = 0.97802556
B = -2.07916483
C = -0.90394038
opts = Odeoptions()
opts.odeset('InitialStep',0.006940920311305178)
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 7
INPUT = INPUT + 1
tspan = [-0.48332446, 14.34676291]
y0_1 = [1.45812793]
y0_2 = [-0.23465422, -3.00879572]
y0_3 = [ 5.20014711,  3.49909892, -2.43879495]
A = -2.52625172
B = -3.07288758
C = 2.71940342
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 8
INPUT = INPUT + 1
opts = Odeoptions()
tspan = [0.5835929 , 11.41241183]
y0_1 = [-6.16053854]
y0_2 = [4.62962776, 6.26526089]
y0_3 = [-6.78523373, -7.61812467, -6.62845498]
A = 4.82456837
B = 4.87407354
C = -0.51263804
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 9
INPUT = INPUT + 1
tspan = [ 0.09661831, 12.72682312]
y0_1 = [-7.86311218]
y0_2 = [2.7171076 , 3.14305633]
y0_3 = [7.4411667 , 5.8554041 , 0.49811558]
A = -1.50925271
B = 4.90586295
C = 2.07144717
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 10
INPUT = INPUT + 1
tspan = [1.16507801, 9.70538009]
y0_1 = [0.16754638]
y0_2 = [-4.10570051,  4.62645345]
y0_3 = [-1.74379021, -5.29778399, -7.08114797]
A = -4.33052396
B = 0.93460987
C = 3.16415395
opts = Odeoptions()
opts.odeset('RelTol',1.927084462770595e-04)
opts.odeset('AbsTol',1.5422113671236964e-06)
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 11
INPUT = INPUT + 1
tspan = [6.06630453, 7.85695321]
y0_1 = [-6.64751259]
y0_2 = [-3.5624895 , -1.92005363]
y0_3 = [ 1.19939346,  5.62223293, -6.47515068]
A = -3.6500774
B = 1.32572362
C = 0.16942695
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 12
INPUT = INPUT + 1
tspan = [9.83453822, 10.33600232]
y0_1 = [-4.2694042]
y0_2 = [7.47046585, 4.87288227]
y0_3 = [-6.23101346, -7.49398607, -2.95299713]
A = -0.16891248
B = 3.4227567
C = 0.33005933
opts = Odeoptions()
opts.odeset('AbsTol',0.00028328386769520237)
opts.odeset('InitialStep',0.03300773645455483)
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 13
INPUT = INPUT + 1 
tspan = [ 3.73942976, 11.92823692]
y0_1 = [0.0485922]
y0_2 = [-2.63842309,  3.56543058]
y0_3 = [-7.42857628,  7.52946529, -1.44846586]
A = 4.60762159
B = 2.90300831
C = 1.71290846
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 14
INPUT = INPUT + 1 
tspan = [0.23664909, 3.47709274]
y0_1 = [5.46724578]
y0_2 = [-7.31855781,  2.52616015]
y0_3 = [ 4.50310439, -6.63021287,  7.36660554]
A = 2.17565639
B = 2.29407597
C = -4.95602666
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 15
INPUT = INPUT + 1 
tspan = [4.26795088, 5.14903329]
y0_1 = [2.25012195]
y0_2 = [7.48015429,  6.35063639]
y0_3 = [4.1022839 , 7.47262896, 0.03892344]
A = 0.53474375
B = 1.07537947
C = 4.95073137
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 16
INPUT = INPUT + 1 
tspan = [ 3.56324789, 12.2673467 ]
y0_1 = [-2.26297594]
y0_2 = [ 2.02273623, -6.45190008]
y0_3 = [-0.10547361,  1.9468299 , -7.04713625]
A = 2.36580799
B = -3.06154501
C = -2.35234947
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 17
INPUT = INPUT + 1 
tspan = [-0.30700883, 12.25575123]
y0_1 = [0.29077201]
y0_2 = [-3.41440131, -6.44950855]
y0_3 = [-1.17725091,  1.92059572, -1.70966996]
A = 2.95104603
B = -0.87055363
C = -3.3465916
opts = Odeoptions()
opts.odeset('RelTol',0.006736304816562963)
opts.odeset('AbsTol',1.7932625474019877e-06)

temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 18
INPUT = INPUT + 1 
tspan = [2.0505969 , 4.07595794]
y0_1 = [7.12694417]
y0_2 = [ 3.56975967, -2.07817288]
y0_3 = [ 7.98461333, -0.59689978, -4.61550035]
A = -3.47937704
B = -3.0336046
C = 2.73797931
opts = Odeoptions()
opts.odeset('RelTol',0.007654990276229565)
opts.odeset('InitialStep',0.0008513865364431812)
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 19
INPUT = INPUT + 1 
tspan = [2.24365274, 8.88873432]
y0_1 = [-2.29474435]
y0_2 = [ 7.41550739, -7.02410812]
y0_3 = [ 7.44402733, -0.55342101, -6.29964052]
A = -1.91085349
B = 4.64813772
C = 2.83354041
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 20
INPUT = INPUT + 1 
tspan = [ 9.46232509, 12.1458127]
y0_1 = [5.7059878]
y0_2 = [-3.07163682,  5.98669745]
y0_3 = [ 2.54796126, -7.00795209,  1.20623878]
A = 3.44085258
B = -4.57285375
C = -3.82489802
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 21
INPUT = INPUT + 1 
tspan = [3.35966743, 6.76027352]
y0_1 = [7.73557472]
y0_2 = [ 2.96211931, 0.53208961]
y0_3 = [5.89125929, 6.10059695,  3.19421586]
A = 3.81728708
B = 2.10322666
C = 3.58740611
opts = Odeoptions()
opts.odeset('AbsTol',0.00016540018580299407)
opts.odeset('InitialStep',0.08911214687593723)
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 22
INPUT = INPUT + 1 
tspan = [ 0.85168538, 11.52706215]
y0_1 = [-5.75267113]
y0_2 = [3.53587685, 3.80035726]
y0_3 = [-7.54188637,  4.09083663, -1.50037043]
A = -3.10994784
B = -2.20786117
C = 3.2113363
opts = Odeoptions()
opts.odeset('RelTol',1.0288230299897688e-6)
opts.odeset('AbsTol',4.84787174455511e-07)
opts.odeset('InitialStep',0.003795204821839059)
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 23
INPUT = INPUT + 1
tspan = [-0.90957528, 10.39029072]
y0_1 = [-5.62286875]
y0_2 = [ 7.96578233, -7.41620033]
y0_3 = [-5.26124403,  5.66296662, -2.2004748 ]
A = -3.17662816
B = -4.87848241
C = 1.08520933
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 24
INPUT = INPUT + 1 
tspan = [3.52514894, 6.89384038]
y0_1 = [0.06673227]
y0_2 = [ 4.06547772, -5.77667385]
y0_3 = [-6.09485915,  2.83317466,  2.84113184]
A = -2.76688306
B = 2.19180466
C = -4.79997074
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 25
INPUT = INPUT + 1 
tspan = [7.57582624, 16.15936357]
y0_1 = [1.49109626]
y0_2 = [-4.56886996, -4.18827044]
y0_3 = [ 7.00094593, -1.76127202, -1.00658487]
A = 1.978158
B = -3.53996799
C = 1.24866477
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 26
INPUT = INPUT + 1 
tspan = [9.29247421, 11.20154662]
y0_1 = [6.42056626]
y0_2 = [ 3.02751907, 7.23128958]
y0_3 = [6.8328569 ,  7.95585027,  5.25184783]
A = 3.6811953
B = 1.09621923
C = 1.53552888
opts = Odeoptions()
opts.odeset('InitialStep',0.0018545051990431387)
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 27
INPUT = INPUT + 1 
tspan = [-4.88825142, 14.02972186]
y0_1 = [6.77769352]
y0_2 = [-0.1580719 , -2.58222842]
y0_3 = [-1.84022282, -4.61121715, -6.75090844]
A = 4.68115446
B = -2.46775989
C = -2.48577196
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 28
INPUT = INPUT + 1 
tspan = [ 7.64008102, 14.19173997]
y0_1 = [3.17532075]
y0_2 = [ 3.03633068, -6.16179753]
y0_3 = [ 6.16866933, -1.81803016, -1.92075269]
A = 4.79308367
B = 2.93275077
C = -2.81334205
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 29
INPUT = INPUT + 1 
tspan = [4.40592519, 6.77211226]
y0_1 = [-1.88405234]
y0_2 = [3.41725858, 2.64320194]
y0_3 = [ 3.32851184,  2.69142253, -5.39137005]
A = -4.89867833
B = 4.9015022
C = -1.20383235
opts = Odeoptions()
opts.odeset('RelTol',1.4648154972812207e-08)
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 30
INPUT = INPUT + 1 
tspan = [-2.29121844, 9.03477936]
y0_1 = [1.7769748]
y0_2 = [-7.99034351,  1.84091706]
y0_3 = [3.97775301, 6.53477137, 6.98213957]
A = 3.48665123
B = -1.78313751
C = -1.70433862
opts = Odeoptions()
opts.odeset('AbsTol',9.284105257511226e-2)
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 31
INPUT = INPUT + 1 
tspan = [3.09505651, 9.03974445]
y0_1 = [3.361784]
y0_2 = [-1.24808648,  3.54409479]
y0_3 = [ 2.11037342, -0.3472534 ,  3.95404699]
A = 3.6168676
B = -2.6285962
C = -0.04739296
opts = Odeoptions()
opts.odeset('RelTol',1.2925459834340685e-06)
opts.odeset('InitialStep',1.062281191880226e-05)
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 32
INPUT = INPUT + 1 
tspan = [0.41345667, 2.24521564]
y0_1 = [-2.8119527]
y0_2 = [2.83055902, 0.41888443]
y0_3 = [7.74890026, -0.01720299,  7.74763048]
A = 4.21513872
B = 0.61825742
C = -2.4283297
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 33
INPUT = INPUT + 1 
tspan = [5.65286451, 6.55329414]
y0_1 = [-6.00439257]
y0_2 = [-7.96559742,  2.72646622]
y0_3 = [ 6.19890213, -2.25338557,  6.09858429]
A = 3.22050966
B = 0.59285462
C = 3.9575617
opts = Odeoptions()
opts.odeset('RelTol',3.464885570872277e-07)
opts.odeset('AbsTol',2.716016891125866e-08)
opts.odeset('InitialStep',0.234944963706823)
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 34
INPUT = INPUT + 1 
tspan = [7.25734918, 11.8754668]
y0_1 = [-1.4258659]
y0_2 = [-5.23934392, -1.76197689]
y0_3 = [-5.65019124, -5.02152462,  7.71309572]
A = -2.10258408
B = -1.58306468
C = 0.67760786
opts = Odeoptions()
opts.odeset('RelTol',0.0007149302009860351)
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 35
INPUT = INPUT + 1 
tspan = [9.69755071, 14.27281151]
y0_1 = [4.91590426]
y0_2 = [-5.53783951, -1.2282395 ]
y0_3 = [-5.35875016, -1.18568651,  0.73626608]
A = -4.68940883
B = 4.87638203
C = 2.92662521
opts = Odeoptions()
temps[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

######################### BAR PLOT

temps10 = 0
temps20 = 0
temps30 = 0
temps40 = 0
temps50 = 0
temps60 = 0
temps70 = 0
temps80 = 0
temps90 = 0
temps100 = 0
temps110 = 0
temps120 = 0
temps130 = 0
temps_plus = 0

for i in range(nbr_Input) :
    for j in range(nbr_test) :
        if temps[i,j] < 0.1 :
            temps10 = temps10 + 1
        elif temps[i,j] < 0.2 :
            temps20 = temps20 +1
        elif temps[i,j] < 0.3 :
            temps30 = temps30 +1
        elif temps[i,j] < 0.4 :
            temps40 = temps40 +1
        elif temps[i,j] < 0.5 :
            temps50 = temps50 +1
        elif temps[i,j] < 0.6 :
            temps60 = temps60 +1
        elif temps[i,j] < 0.7 :
            temps70 = temps70 +1
        elif temps[i,j] < 0.8 :
            temps80 = temps80 +1
        elif temps[i,j] < 0.9 :
            temps90 = temps90 +1
        elif temps[i,j] < 1 :
            temps100 = temps100 +1
        elif temps[i,j] < 1.1 :
            temps110 = temps110 +1
        elif temps[i,j] < 1.2 :
            temps120 = temps120 +1
        elif temps[i,j] < 1.3 :
            temps130 = temps130 +1
        else :
            temps_plus = temps_plus +1
            
print([temps10, temps20, temps30, temps40, temps50, temps60, temps70, temps80, temps90, temps100, temps110, temps120, temps130, temps_plus])
#Make barplot

height = [temps30, temps40, temps50, temps60, temps70, temps80, temps90, temps100, temps110, temps120, temps130, temps_plus]
bars = ('< 30%', '< 40%', '< 50%', '< 60%', '< 70%', '< 80%', '< 90%', '< 100%', '< 110%', '< 120%', '< 130%', '> 130%')
y_pos = np.arange(len(bars))

# Create bars
plt.bar(y_pos, height)

# Create names on the x-axis
plt.xticks(y_pos, bars)

# Show graphic
plt.show()