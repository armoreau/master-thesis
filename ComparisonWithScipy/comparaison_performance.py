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
from ode45 import ode45

######################################### function test
### First order fction
### trigonometrique
def f1(t,y) :
    dydt = np.cos(t)
    return dydt

def solf1(t,t0,y0):
    return np.sin(t) + y0 -np.sin(t0)

def f2(t,y) :
    dydt = np.sin(t)
    return dydt

def solf2(t,t0,y0):
    return -np.cos(t) + y0 + np.cos(t0)

def f3(t,y) :
    dydt = np.sin(t) - np.cos(t)
    return dydt

def solf3(t,t0,y0):
    return - np.sin(t) -np.cos(t) + np.cos(t0) + np.sin(t0) + y0

### exponential
def f4(t,y) :
    dydt = np.exp(t)
    return dydt

def solf4(t,t0,y0):
    return np.exp(t) + y0 - np.exp(t0)

def f5(t,y) :
    dydt = y
    return dydt

def solf5(t,t0,y0):
    return y0*np.exp(t-t0)

def f6(t,y) :
    dydt = [-0.5*y]
    return dydt

def solf6(t,t0,y0):
    return y0*np.exp(1)**(0.5*t0-0.5*t)

def f7(t,y) :
    dydt = 0.5*y
    return dydt

def solf7(t,t0,y0):
    return y0*np.exp(1)**(-0.5*t0+0.5*t)

### polynomial
def f8(t,y) :
    dydt = 2*t + np.exp(t)
    return dydt

def solf8(t,t0,y0):
    return t**2 - t0**2 + y0 + np.exp(t) - np.exp(t0)

def f9(t,y) :
    dydt = t**2 + 2*t + np.exp(-t)
    return dydt

def solf9(t,t0,y0):
    return - (t0**3)/3 - t0**2 + y0 + (t**3)/3 + t**2 - np.exp(-t) + np.exp(-t0)

def f10(t,y) :
    dydt = 5*t**5 - t**2 + 2*t 
    return dydt

def solf10(t,t0,y0):
    return (1/6)*(-5*t0**6 + 2*t0**3 - 6*t0**2 + 6*y0 + (t**2) * (5*t**4 - 2*t + 6))

### other
def f11(t,y) :
    dydt = y + t
    return dydt

def solf11(t,t0,y0):
    return (t0 + y0 + 1)*np.exp(t-t0) - t - 1

def f12(t,y) :
    dydt = np.sin(t) + t
    return dydt

def solf12(t,t0,y0):
    return - (t0**2)/2 + np.cos(t0) + y0 + (t**2)/2 - np.cos(t)

def f13(t,y) :
    dydt = np.sin(t) + t**2
    return dydt

def solf13(t,t0,y0):
    return -(t0**3)/3 + np.cos(t0) + y0 + (t**3)/3 - np.cos(t)

def f14(t,y) :
    dydt = 2*y+t**2
    return dydt

def solf14(t,t0,y0):
    return (1/4)*((2*t0**2 + 2*t0 + 4*y0 +1)*np.exp(2*t - 2*t0) - 2*t**2 -2*t - 1) 

def f15(t,y) :
    dydt = np.sin(t) + 1/t
    return dydt

def solf15(t,t0,y0):
    return -np.log(t0) + np.cos(t0) + y0 + np.log(t) - np.cos(t)

#### second order or more fction

def f16(t,y) :
    dydt = np.zeros(2)
    dydt[0] = t
    dydt[1] = y[1]
    return dydt

def solf16(t,t0,y0):
    sol = np.zeros(2)
    sol[0] = y0[0] - (t0**2)/2 + (t**2)/2
    sol[1] = (y0[1]*np.exp(-t0))*np.exp(t)
    return sol

def f17(t,y) :
    dydt = np.zeros(2)
    dydt[0] = y[1]
    dydt[1] = y[0]
    return dydt

def solf17(t,t0,y0): #t0 = 0
    sol = np.zeros(2)
    sol[0] = (1/2) * np.exp(-t)*(y0[0]*(np.exp(2*t)+1) + y0[1]*(np.exp(2*t)-1))
    sol[1] = (1/2) * np.exp(-t)*(y0[0]*(np.exp(2*t)-1) + y0[1]*(np.exp(2*t)+1))
    return sol

def f18(t,y) :
    dydt = np.zeros(2)
    dydt[0] = 2*y[0]
    dydt[1] = 1/(y[0]+100)
    return dydt

def solf18(t,t0,y0): #t0 = 0
    sol = np.zeros(2)
    sol[0] = y0[0]*np.exp(2*t)
    sol[1] = (1/200)*(-np.log(y0[0]*np.exp(2*t) +100) + np.log(y0[0] +100) +200*y0[1] +2*t)
    return sol
#### additionnal argument function

def f19(t,y,A,B,C) :
    dydt = (A+B)*y/C
    return dydt

def solf19(t,t0,y0,A,B,C):
    return y0*np.exp(-((t0-t)*(A+B))/C)

def f20(t,y,A,B,C) :
    dydt = np.sin(t*(A+B)/C)
    return dydt

def solf20(t,t0,y0,A,B,C):
    return (C*np.cos((t0*(A+B))/C) + y0*(A+B) - C*np.cos((t*(A+B))/C))/(A+B)

def compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C) :
    
    precision_mean = np.zeros([20])
    j=0
    t0 = tspan[0]
    
    y0 = np.array(y0_1)
    arg = None
    
    scipy_sol = solve_ivp(f1,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f1, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf1(tspan[i],t0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    scipy_sol = solve_ivp(f2,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f2, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf2(tspan[i],t0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    scipy_sol = solve_ivp(f3,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f3, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf3(tspan[i],t0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    scipy_sol = solve_ivp(f4,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f4, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf4(tspan[i],t0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    scipy_sol = solve_ivp(f5,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f5, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf5(tspan[i],t0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    scipy_sol = solve_ivp(f6,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f6, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf6(tspan[i],t0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    scipy_sol = solve_ivp(f7,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f7, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf7(tspan[i],t0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    scipy_sol = solve_ivp(f8,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f8, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf8(tspan[i],t0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    scipy_sol = solve_ivp(f9,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f9, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf9(tspan[i],t0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    scipy_sol = solve_ivp(f10,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f10, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf10(tspan[i],t0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    scipy_sol = solve_ivp(f11,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f11, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf11(tspan[i],t0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    scipy_sol = solve_ivp(f12,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f12, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf12(tspan[i],t0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    scipy_sol = solve_ivp(f13,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f13, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf13(tspan[i],t0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    scipy_sol = solve_ivp(f14,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f14, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf14(tspan[i],t0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    scipy_sol = solve_ivp(f15,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f15, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf15(tspan[i],t0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    ### y0
    y0 = np.array(y0_2)
    
    scipy_sol = solve_ivp(f16,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f16, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf16(tspan[i],t0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    #### t0 must be 0
    tspan2 = np.linspace(0,12,20)
    
    scipy_sol = solve_ivp(f17,[tspan2[0],tspan2[-1]],y0,t_eval=tspan2,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f17, tspan2, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan2)])
    for i in range(len(tspan2)) :
        true_sol[:,i] = solf17(tspan2[i],0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    scipy_sol = solve_ivp(f18,[tspan2[0],tspan2[-1]],y0,t_eval=tspan2,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f18, tspan2, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan2)])
    for i in range(len(tspan2)) :
        true_sol[:,i] = solf18(tspan2[i],0,y0)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    ### y0 
    y0 = np.array(y0_1)
    arg = A,B,C
    
    scipy_sol = solve_ivp(f19,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f19, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf19(tspan[i],t0,y0,A,B,C)
    ########### Compute erreur    
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1
    
    scipy_sol = solve_ivp(f20,[tspan[0],tspan[-1]],y0,t_eval=tspan,args=arg,rtol=opts.RelTol,atol=opts.AbsTol,first_step=opts.InitialStep)
    sol = ode45(f20, tspan, y0, opts, arg)
    ########### True sol
    true_sol = np.zeros([len(y0),len(tspan)])
    for i in range(len(tspan)) :
        true_sol[:,i] = solf20(tspan[i],t0,y0,A,B,C)
    ########### Compute erreur
    abs_error = np.abs(true_sol-sol.y)
    scipy_abs_error = np.abs(true_sol-scipy_sol.y)
    precision_mean[j] = np.mean(scipy_abs_error)/np.mean(abs_error)
    j = j+1

    return precision_mean

######################################### perform test with random input
nbr_test = 20
nbr_Input = 35
precision_mean = np.zeros([nbr_Input, nbr_test])

#### INPUT 1
INPUT = 0
tspan = [ 0.74414985,  3.5823718 ,  6.42059375,  9.2588157 , 12.09703764,
       14.93525959, 17.77348154, 20.61170349, 23.44992544, 26.28814739,
       29.12636934, 31.96459129]
y0_1 = [4.53897499]
y0_2 = [1.20051259, 1.17666892]
y0_3 = [ 0.45052893,  7.03086018, -1.50526611]
A = -2.83009697
B = -4.97673659
C = -3.50342477
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 2
INPUT = INPUT + 1
tspan = [ 1.69850555,  4.63711663,  7.57572772, 10.5143388 , 13.45294988,
       16.39156097, 19.33017205, 22.26878314, 25.20739422, 28.1460053 ,
       31.08461639, 34.02322747, 36.96183856, 39.90044964, 42.83906072,
       45.77767181, 48.71628289]
y0_1 = [-6.38784938]
y0_2 = [1.60628946, 1.18921585]
y0_3 = [-7.09018137,  4.81021195, -7.82701511]
A = 2.25480317
B = -0.94979655
C = 0.1968337
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)


#### INPUT 3
INPUT = INPUT + 1
tspan =  [ 8.96715752,  9.36035427,  9.75355102, 10.14674776, 10.53994451,
       10.93314125, 11.326338  , 11.71953475, 12.11273149, 12.50592824,
       12.89912498, 13.29232173, 13.68551848, 14.07871522, 14.47191197,
       14.86510871, 15.25830546, 15.6515022 , 16.04469895]
y0_1 = [-1.94429498]
y0_2 = [4.96831322, -7.15530037]
y0_3 = [-2.33703264, -2.00682155, -7.14666703]
A = 4.75708095
B = 3.95418388
C = -4.72783249
opts = Odeoptions()
opts.odeset('RelTol',1.3818760883438436e-04)
opts.odeset('AbsTol',2.0635002062377786e-05)
opts.odeset('InitialStep',0.23252681078089357)
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 4
INPUT = INPUT + 1
tspan = [ 8.39586318,  8.91985661,  9.01857804,  9.08393703,  9.13868348,
       10.11329884, 10.40673529, 10.51920706, 10.99127963, 11.14986855,
       11.21056473, 11.35003671, 11.8351665 , 12.08577212, 12.74535624,
       12.86117488]
y0_1 = [-3.76904574]
y0_2 = [5.24019054,  4.65499285]
y0_3 = [-6.51846199,  3.11856858,  1.27300755]
A = 1.30937634
B = -1.17558159
C = -1.51165755
opts = Odeoptions()
opts.odeset('AbsTol',1.4482218440447706e-09)
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 5
INPUT = INPUT + 1
tspan = [ 8.26274166, 13.2839027 , 18.30506374, 23.32622478, 28.34738582,
       33.36854686, 38.3897079 ]
y0_1 = [-0.01349065]
y0_2 = [0.76687766, 3.93978074]
y0_3 = [2.98288271, 3.48166087, 7.2604428 ]
A = 6.61349864e-01
B = -3.65342548e-03
C = -4.36285600e+00
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 6
INPUT = INPUT + 1
tspan = [ 5.22682674,  6.27095958,  7.31509243,  8.35922527,  9.40335811,
       10.44749095, 11.4916238 , 12.53575664, 13.57988948, 14.62402232,
       15.66815517, 16.71228801, 17.75642085, 18.80055369, 19.84468654,
       20.88881938, 21.93295222, 22.97708506, 24.0212179 ]
y0_1 = [-2.51324011]
y0_2 = [6.18446334,  7.34621481]
y0_3 = [-4.52846858,  2.46373529,  2.13908458]
A = 0.97802556
B = -2.07916483
C = -0.90394038
opts = Odeoptions()
opts.odeset('InitialStep',0.006940920311305178)
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 7
INPUT = INPUT + 1
tspan = [ 9.39678878, 11.65462621, 13.91246364, 16.17030107, 18.4281385 ,
       20.68597593, 22.94381335, 25.20165078, 27.45948821]
y0_1 = [1.45812793]
y0_2 = [0.23465422, 3.00879572]
y0_3 = [ 5.20014711,  3.49909892, -2.43879495]
A = -2.52625172
B = -3.07288758
C = 2.71940342
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 8
INPUT = INPUT + 1
opts = Odeoptions()
tspan = [ 9.81653292, 10.17379192, 11.06642541, 11.27627208, 12.24407322,
       12.33738357, 12.42494572, 13.33218639, 14.17357666, 14.69740422,
       14.78054436, 15.32640978]
y0_1 = [-6.16053854]
y0_2 = [4.62962776, 6.26526089]
y0_3 = [-6.78523373, -7.61812467, -6.62845498]
A = 4.82456837
B = 4.87407354
C = -0.51263804
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 9
INPUT = INPUT + 1
tspan = [ 9.02577645, 10.87750684, 12.72923722, 14.5809676 , 16.43269798,
       18.28442836, 20.13615875]
y0_1 = [-7.86311218]
y0_2 = [2.7171076 , 3.14305633]
y0_3 = [7.4411667 , 5.8554041 , 0.49811558]
A = -1.50925271
B = 4.90586295
C = 2.07144717
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 10
INPUT = INPUT + 1
tspan = [ 5.08711862,  5.81515153,  6.54318445,  7.27121737,  7.99925029,
        8.72728321,  9.45531612, 10.18334904]
y0_1 = [0.16754638]
y0_2 = [4.10570051,  4.62645345]
y0_3 = [-1.74379021, -5.29778399, -7.08114797]
A = -4.33052396
B = 0.93460987
C = 3.16415395
opts = Odeoptions()
opts.odeset('RelTol',1.927084462770595e-04)
opts.odeset('AbsTol',1.5422113671236964e-06)
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 11
INPUT = INPUT + 1
tspan = [ 8.4727269 , 11.93021544, 15.38770398, 18.84519253, 22.30268107,
       25.76016961, 29.21765815, 32.67514669, 36.13263523, 39.59012377,
       43.04761231, 46.50510085, 49.9625894 ]
y0_1 = [-6.64751259]
y0_2 = [3.5624895 , 1.92005363]
y0_3 = [ 1.19939346,  5.62223293, -6.47515068]
A = -3.6500774
B = 1.32572362
C = 0.16942695
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 12
INPUT = INPUT + 1
tspan = [ 3.17345345,  5.40886796,  7.64428247,  9.87969698, 12.11511149,
       14.35052601, 16.58594052, 18.82135503, 21.05676954, 23.29218405,
       25.52759856, 27.76301307]
y0_1 = [-4.2694042]
y0_2 = [7.47046585, 4.87288227]
y0_3 = [-6.23101346, -7.49398607, -2.95299713]
A = -0.16891248
B = 3.4227567
C = 0.33005933
opts = Odeoptions()
opts.odeset('AbsTol',0.00028328386769520237)
opts.odeset('InitialStep',0.03300773645455483)
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 13
INPUT = INPUT + 1 
tspan = [ 9.7453995 ,  9.95090191, 10.67001893, 10.8323155 , 11.20609249,
       11.74565997, 12.47967583, 13.22721628, 14.00644958, 14.03028747,
       14.58681226, 15.10279257, 15.34011124, 16.06370832, 16.50484584,
       17.44923163, 17.6368765 ]
y0_1 = [0.0485922]
y0_2 = [2.63842309,  3.56543058]
y0_3 = [-7.42857628,  7.52946529, -1.44846586]
A = 4.60762159
B = 2.90300831
C = 1.71290846
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 14
INPUT = INPUT + 1 
tspan = [ 2.7622757 , 10.43626758, 18.11025946, 25.78425133, 33.45824321,
       41.13223509]
y0_1 = [5.46724578]
y0_2 = [7.31855781,  2.52616015]
y0_3 = [ 4.50310439, -6.63021287,  7.36660554]
A = 2.17565639
B = 2.29407597
C = -4.95602666
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 15
INPUT = INPUT + 1 
tspan = [ 6.9850128 ,  7.97523936,  8.52096614,  8.88686077,  9.70013126,
        9.81796367,  9.97957887, 10.48182203, 11.04127875, 11.64466402,
       11.67150033, 12.52149521, 13.24222898, 13.58032784, 13.83109881,
       14.37548757]
y0_1 = [2.25012195]
y0_2 = [7.48015429,  6.35063639]
y0_3 = [4.1022839 , 7.47262896, 0.03892344]
A = 0.53474375
B = 1.07537947
C = 4.95073137
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 16
INPUT = INPUT + 1 
tspan = [ 9.13273478,  9.49347325,  9.77797685,  9.91093543, 10.52953047,
       11.22776315, 12.02425432, 12.45007967, 12.45928544, 13.06569571,
       13.43964931, 13.50069978, 14.03579572]
y0_1 = [-2.26297594]
y0_2 = [ 2.02273623, 6.45190008]
y0_3 = [-0.10547361,  1.9468299 , -7.04713625]
A = 2.36580799
B = -3.06154501
C = -2.35234947
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 17
INPUT = INPUT + 1 
tspan = [ 5.40360036,  5.9946027 ,  6.11202109,  7.00658649,  7.36171297,
        7.42131194,  8.29171667,  9.19710614,  9.65094862, 10.18197877,
       11.17077419, 11.21934425, 11.29696039, 11.63429422, 12.49039877,
       12.84653281]
y0_1 = [0.29077201]
y0_2 = [3.41440131, 6.44950855]
y0_3 = [-1.17725091,  1.92059572, -1.70966996]
A = 2.95104603
B = -0.87055363
C = -3.3465916
opts = Odeoptions()
opts.odeset('RelTol',0.006736304816562963)
opts.odeset('AbsTol',1.7932625474019877e-06)

precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 18
INPUT = INPUT + 1 
tspan = [ 4.87234942,  5.00440906,  5.64415681,  5.76514867,  6.45246345,
        7.14854372,  7.96956847,  8.77365716,  9.12524841,  9.82213236,
       10.03750046, 10.41505654, 10.5758234 ]
y0_1 = [7.12694417]
y0_2 = [ 3.56975967, 2.07817288]
y0_3 = [ 7.98461333, -0.59689978, -4.61550035]
A = -3.47937704
B = -3.0336046
C = 2.73797931
opts = Odeoptions()
opts.odeset('RelTol',0.007654990276229565)
opts.odeset('InitialStep',0.0008513865364431812)
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 19
INPUT = INPUT + 1 
tspan = [ 3.08131047,  3.14105644,  3.55677216,  3.92552722,  4.81645502,
        5.51844234,  6.38882331,  7.34239687,  7.98908828,  8.43820397,
        9.28342975,  9.82522026,  9.89049811, 10.86556611, 11.71479577,
       12.27775326, 13.00095614, 13.72914177]
y0_1 = [-2.29474435]
y0_2 = [ 7.41550739, 7.02410812]
y0_3 = [ 7.44402733, -0.55342101, -6.29964052]
A = -1.91085349
B = 4.64813772
C = 2.83354041
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 20
INPUT = INPUT + 1 
tspan = [ 2.91402466,  3.41863076,  4.27858316,  4.61503076,  4.77487416,
        5.24610874,  5.72983883,  6.44321663,  6.64099038,  7.42071527,
        7.80825897,  8.68595965,  9.10784503,  9.58161657, 10.00677973]
y0_1 = [5.7059878]
y0_2 = [3.07163682,  5.98669745]
y0_3 = [ 2.54796126, -7.00795209,  1.20623878]
A = 3.44085258
B = -4.57285375
C = -3.82489802
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 21
INPUT = INPUT + 1 
tspan = [ 4.26233857,  5.75510873,  7.24787888,  8.74064903, 10.23341918,
       11.72618933, 13.21895949, 14.71172964, 16.20449979, 17.69726994,
       19.19004009, 20.68281025]
y0_1 = [7.73557472]
y0_2 = [ 2.96211931, 0.53208961]
y0_3 = [5.89125929, 6.10059695,  3.19421586]
A = 3.81728708
B = 2.10322666
C = 3.58740611
opts = Odeoptions()
opts.odeset('AbsTol',0.00016540018580299407)
opts.odeset('InitialStep',0.08911214687593723)
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 22
INPUT = INPUT + 1 
tspan = [ 8.557499  , 10.76548469, 12.97347039, 15.18145609, 17.38944178,
       19.59742748, 21.80541318, 24.01339888, 26.22138457, 28.42937027,
       30.63735597, 32.84534167, 35.05332736, 37.26131306, 39.46929876,
       41.67728446, 43.88527015, 46.09325585]
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
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 23
INPUT = INPUT + 1
tspan = [ 1.61245558,  2.20996442,  2.51457378,  3.08621419,  3.73241815,
        4.33745399,  5.32112185,  6.21885974,  6.68054178,  7.34299055,
        8.15589473,  8.34430679,  8.4820539 ,  8.73776518,  9.48109092,
       10.03120337, 10.88030284, 11.14367106, 11.75764566]
y0_1 = [-5.62286875]
y0_2 = [ 7.96578233, 7.41620033]
y0_3 = [-5.26124403,  5.66296662, -2.2004748 ]
A = -3.17662816
B = -4.87848241
C = 1.08520933
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 24
INPUT = INPUT + 1 
tspan = [ 4.08250894,  5.63894718,  7.19538542,  8.75182365, 10.30826189,
       11.86470013, 13.42113837, 14.9775766 , 16.53401484, 18.09045308,
       19.64689132, 21.20332956, 22.75976779, 24.31620603, 25.87264427,
       27.42908251, 28.98552074, 30.54195898]
y0_1 = [0.06673227]
y0_2 = [ 4.06547772, -5.77667385]
y0_3 = [-6.09485915,  2.83317466,  2.84113184]
A = -2.76688306
B = 2.19180466
C = -4.79997074
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 25
INPUT = INPUT + 1 
tspan = [ 5.08159835,  5.86084285,  6.1134522 ,  7.0071192 ,  7.53890476,
        8.13023418,  8.14911493,  8.96205468,  9.54978561,  9.66938852,
       10.43986065, 11.04564893]
y0_1 = [1.49109626]
y0_2 = [4.56886996, 4.18827044]
y0_3 = [ 7.00094593, -1.76127202, -1.00658487]
A = 1.978158
B = -3.53996799
C = 1.24866477
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 26
INPUT = INPUT + 1 
tspan = [ 7.73226125,  9.16100248, 10.5897437 , 12.01848493, 13.44722616,
       14.87596739, 16.30470861, 17.73344984, 19.16219107, 20.59093229,
       22.01967352, 23.44841475, 24.87715598, 26.3058972 , 27.73463843,
       29.16337966]
y0_1 = [6.42056626]
y0_2 = [ 3.02751907, 7.23128958]
y0_3 = [6.8328569 ,  7.95585027,  5.25184783]
A = 3.6811953
B = 1.09621923
C = 1.53552888
opts = Odeoptions()
opts.odeset('InitialStep',0.0018545051990431387)
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 27
INPUT = INPUT + 1 
tspan = [ 9.19397481,  9.72971727, 10.36474027, 11.18196793, 11.76297821,
       11.84968782, 12.17164847, 12.22664662, 13.20482502, 13.9492462 ,
       14.54591319, 14.72382504, 15.15886832, 15.57880606, 16.35276041,
       16.59579624]
y0_1 = [6.77769352]
y0_2 = [0.1580719 , 2.58222842]
y0_3 = [-1.84022282, -4.61121715, -6.75090844]
A = 4.68115446
B = -2.46775989
C = -2.48577196
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 28
INPUT = INPUT + 1 
tspan = [ 8.50675232,  8.92145336,  9.41487621, 10.29223962, 10.84873283,
       10.94359994, 10.99164235, 11.82099971, 12.06582897, 12.49390149,
       12.5206252 , 12.69962281, 13.32476178, 14.25650939, 14.69033935,
       14.94685873]
y0_1 = [3.17532075]
y0_2 = [ 3.03633068, 6.16179753]
y0_3 = [ 6.16866933, -1.81803016, -1.92075269]
A = 4.79308367
B = 2.93275077
C = -2.81334205
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 29
INPUT = INPUT + 1 
tspan = [ 8.14471516, 10.37006175, 12.59540834, 14.82075493, 17.04610152,
       19.27144811, 21.49679471, 23.7221413 , 25.94748789, 28.17283448,
       30.39818107, 32.62352766]
y0_1 = [-1.88405234]
y0_2 = [3.41725858, 2.64320194]
y0_3 = [ 3.32851184,  2.69142253, -5.39137005]
A = -4.89867833
B = 4.9015022
C = -1.20383235
opts = Odeoptions()
opts.odeset('RelTol',1.4648154972812207e-08)
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 30
INPUT = INPUT + 1 
tspan = [ 8.01517385,  8.97569477,  9.42050561,  9.92138412, 10.06690694,
       10.45460128, 10.71547724, 11.54870089, 12.3332083 , 12.64375186,
       13.01052409, 13.01925285, 13.0749965 , 13.36858153, 14.14815074,
       14.47086066]
y0_1 = [1.7769748]
y0_2 = [7.99034351,  1.84091706]
y0_3 = [3.97775301, 6.53477137, 6.98213957]
A = 3.48665123
B = -1.78313751
C = -1.70433862
opts = Odeoptions()
opts.odeset('AbsTol',9.284105257511226e-2)
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 31
INPUT = INPUT + 1 
tspan = [6.4588985 , 6.58767103, 6.71644355, 6.84521608, 6.9739886 ,
       7.10276113, 7.23153365, 7.36030618, 7.4890787 , 7.61785123,
       7.74662375, 7.87539628, 8.0041688 , 8.13294133, 8.26171385,
       8.39048638, 8.5192589 ]
y0_1 = [3.361784]
y0_2 = [1.24808648,  3.54409479]
y0_3 = [ 2.11037342, -0.3472534 ,  3.95404699]
A = 3.6168676
B = -2.6285962
C = -0.04739296
opts = Odeoptions()
opts.odeset('RelTol',1.2925459834340685e-02)
opts.odeset('InitialStep',1.062281191880226e-05)
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 32
INPUT = INPUT + 1 
tspan = [ 3.26822329,  3.7486273 ,  4.22903131,  4.70943533,  5.18983934,
        5.67024336,  6.15064737,  6.63105138,  7.1114554 ,  7.59185941,
        8.07226342,  8.55266744,  9.03307145,  9.51347547,  9.99387948,
       10.47428349, 10.95468751, 11.43509152]
y0_1 = [-2.8119527]
y0_2 = [2.83055902, 0.41888443]
y0_3 = [7.74890026, -0.01720299,  7.74763048]
A = 4.21513872
B = 0.61825742
C = -2.4283297
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 33
INPUT = INPUT + 1 
tspan = [2.71797756, 3.11835641, 3.71726714, 4.68239104, 5.49341035,
       6.42120491, 6.84890184, 6.89239923, 7.68396108, 8.58018338,
       9.52656294, 9.57933574]
y0_1 = [-6.00439257]
y0_2 = [7.96559742,  2.72646622]
y0_3 = [ 6.19890213, -2.25338557,  6.09858429]
A = 3.22050966
B = 0.59285462
C = 3.9575617
opts = Odeoptions()
opts.odeset('RelTol',3.464885570872277e-02)
opts.odeset('AbsTol',2.716016891125866e-02)
opts.odeset('InitialStep',0.234944963706823)
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 34
INPUT = INPUT + 1 
tspan = [ 3.9934111 ,  4.85325873,  5.26668727,  6.09458822,  6.57793317,
        6.68838707,  6.95377958,  7.07125618,  7.88009837,  8.8145981 ,
        8.99520036,  9.26528354,  9.93840577, 10.50059607, 10.67100908]
y0_1 = [-1.4258659]
y0_2 = [5.23934392, 1.76197689]
y0_3 = [-5.65019124, -5.02152462,  7.71309572]
A = -2.10258408
B = -1.58306468
C = 0.67760786
opts = Odeoptions()
opts.odeset('RelTol',0.0007149302009860351)
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

#### INPUT 35
INPUT = INPUT + 1 
tspan = [1.94698349, 2.43500423, 3.04980331, 3.93601644, 4.14004799,
       4.57325309, 5.53196007, 5.91691221, 6.29274032, 6.92900246,
       7.24446405, 8.19971252, 8.95330826]
y0_1 = [4.91590426]
y0_2 = [5.53783951, 1.2282395]
y0_3 = [-5.35875016, -1.18568651,  0.73626608]
A = -4.68940883
B = 4.87638203
C = 2.92662521
opts = Odeoptions()
precision_mean[INPUT,:] = compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C)

######################### BAR PLOT

precision001 = 0
precision01 = 0
precision1 = 0
precision10 = 0
precision100 = 0
precision1000 = 0
precision10000 = 0
precision100000 = 0
precision_moins = 0


for i in range(nbr_Input) :
    for j in range(nbr_test) :
        if precision_mean[i,j] > 100000 :
            precision100000 = precision100000 + 1
        elif precision_mean[i,j] > 10000 :
            precision10000 = precision10000 +1
        elif precision_mean[i,j] > 1000 :
            precision1000 = precision1000 +1
        elif precision_mean[i,j] > 100 :
            precision100 = precision100 +1
        elif precision_mean[i,j] > 10 :
            precision10 = precision10 +1
        elif precision_mean[i,j] > 1 :
            precision1 = precision1 +1
        elif precision_mean[i,j] > 0.1 :
            precision01 = precision01 +1
        elif precision_mean[i,j] > 0.01 :
            precision001 = precision001 +1
        else :
            precision_moins = precision_moins +1
            
print([precision100000, precision10000, precision1000, precision100, precision10, precision1, precision01, precision001, precision_moins])
#Make barplot

height = [precision100000, precision10000, precision1000, precision100, precision10, precision1, precision01, precision001, precision_moins]
bars = ('> 100000', '> 10000', '> 1000', '> 100', '> 10', '> 1', '> 0.1', '>0.01', '0.001')
y_pos = np.arange(len(bars))

# Create bars
plt.bar(y_pos, height)

# Create names on the x-axis
plt.xticks(y_pos, bars)

# Show graphic
plt.show()