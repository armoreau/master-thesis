from scipy.integrate import solve_ivp
from ode45 import ode45
from Options import Options
import numpy as np
import matplotlib.pyplot as plt
import time

################ Function test
def fun1(t, y):
    dydt = [-0.5*y]
    return dydt

def fun2(t,y) :
    dydt = np.cos(t)
    return dydt

def fun3(t, y):
    dydt = np.array([y[1]+y[2], -y[0]+2*y[1]+y[2], y[0]+y[2]])
    return dydt

def fun4(t,y) :
    dydt = np.zeros(2)
    dydt[0] = y[1]
    dydt[1] = (1-y[0]**2)*y[1]-y[0]
    return dydt

def fun5(t,y,A,B) :
    dydt = np.array([y[1], (A/B)*t* y[0]])
    return dydt

def fun6(t,y) :
    dydt = -2*y + 2*np.cos(t)*np.sin(2*t)
    return dydt

def fun7(t,y) :
    mu = 1 / 82.45
    mustar = 1 - mu
    r13 = ((y[0] + mu)**2 + y[1]**2) ** 1.5
    r23 = ((y[0] - mustar)**2 + y[1]**2) ** 1.5
    dydt = np.array([ y[2],y[3],2*y[3] + y[0] - mustar*((y[0]+mu)/r13) - mu*((y[0]-mustar)/r23),-2*y[2] + y[1] - mustar*(y[1]/r13) - mu*(y[1]/r23)])
    return dydt

def events_fun7(t,y):
    y0 = np.array([1.2, 0, 0, -1.04935750])
    dDSQdt = 2 * np.dot((y[0:2]-y0[0:2]),y[2:4])
    value = np.array([dDSQdt, dDSQdt])
    isterminal = np.array([1, 0])
    direction = np.array([1, -1])
    return [value,isterminal,direction]

def fun8(t,y) :
    m1 = 0.1
    m2 = 0.1
    L = 1
    g = 9.81
    dydt = np.array([y[1], m2*L*(y[5]**2)*np.cos(y[4]),y[3],m2*L*(y[5]**2)*np.sin(y[4])-(m1+m2)*g,y[5],-g*L*np.cos(y[4])])
    return dydt

def mass_fun8(t,y) :
    m1 = 0.1
    m2 = 0.1
    L = 1
    M = np.zeros([6,6])
    M[0,0] = 1
    M[1,1] = m1 + m2
    M[1,5] = -m2*L*np.sin(y[4])
    M[2,2] = 1
    M[3,3] = m1 + m2
    M[3,5] = m2*L*np.cos(y[4])
    M[4,4] = 1
    M[5,1] = -L*np.sin(y[4])
    M[5,3] = L*np.cos(y[4])
    M[5,5] = L**2
    return M

def fun9(t,y) :
    dydt = np.array([y[1],-5*y[1]+4*y[0]+np.sin(10*t)])
    return dydt

######################################### Test 1

tspan = np.linspace(0,10,25)
y0 = np.array([2])

sol = ode45(fun1, tspan, y0)

np.savetxt('Test1_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test1_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 2

tspan = np.linspace(0,10,25)
y0 = np.array([0])

sol = ode45(fun2, tspan, y0)

np.savetxt('Test2_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test2_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 3

tspan = np.linspace(0,2,25)
y0 = np.array([2, 4, 6])

sol = ode45(fun3, tspan, y0)

np.savetxt('Test3_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test3_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 4

tspan = np.linspace(0,20,25)
y0 = np.array([2, 0])

sol = ode45(fun4, tspan, y0)

np.savetxt('Test4_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test4_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 5

A=1
B=2
arg = A,B
tspan = np.linspace(0,5,25)
y0 = np.array([0, 0.01])

sol = ode45(fun5, tspan, y0, varargin = arg)

np.savetxt('Test5_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test5_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 6

tspan = np.linspace(0,3,25)
y0 = np.array([-3,-2,-1,0,1,2,3])

sol = ode45(fun6, tspan, y0)

np.savetxt('Test6_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test6_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 7
tspan = np.linspace(0,10,25)
y0 = np.array([2])

myOptions = Options() #Create an option with default value.
myOptions.odeset('Refine',10)
myOptions.odeset('NormControl',True)
myOptions.odeset('MaxStep',1)
myOptions.odeset('InitialStep',np.array([0.1]))

sol = ode45(fun1,tspan,y0, myOptions)

np.savetxt('Test7_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test7_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 8
tspan = np.linspace(0,5,25)
y0 = np.array([0, 0.01])
A=1
B=2
arg = A,B

myOptions = Options() #Create an option with default value.
myOptions.odeset('Refine',1)
myOptions.odeset('NormControl',False)
myOptions.odeset('MaxStep',np.array([0.5]))
myOptions.odeset('InitialStep',np.array([0.1]))

sol = ode45(fun5,tspan,y0, myOptions, arg)

np.savetxt('Test8_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test8_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 9

tspan = np.linspace(0,3,25)
y0 = np.array([0, 0])

sol = ode45(fun9,tspan,y0)

np.savetxt('Test9_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test9_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 10

tspan = np.linspace(0,4,25)
y0 = [0, 4, 2, 20, -np.pi/2, 2]
options = Options()
options.odeset('Mass',mass_fun8)

sol = ode45(fun8,tspan,y0,options)

np.savetxt('Test10_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test10_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 11

tspan = np.linspace(10,0,25)
y0 = np.array([2])

sol = ode45(fun1, tspan, y0)

np.savetxt('Test11_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test11_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 12

tspan = np.linspace(10,0,25)
y0 = np.array([0])

sol = ode45(fun2, tspan, y0)

np.savetxt('Test12_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test12_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 13

tspan = np.linspace(2,0,25)
y0 = np.array([2, 4, 6])

sol = ode45(fun3, tspan, y0)

np.savetxt('Test13_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test13_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 14

tspan = np.linspace(20,0,25)
y0 = np.array([2, 0])

sol = ode45(fun4, tspan, y0)

np.savetxt('Test14_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test14_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 15

A=1
B=2
arg = A,B
opts = Options()
tspan = np.linspace(5,0,25)
y0 = np.array([0, 0.01])

sol = ode45(fun5, tspan, y0,opts,arg)

np.savetxt('Test15_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test15_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 16

tspan = np.linspace(3,0,25)
y0 = np.array([-3,-2,-1,0,1,2,3])

sol = ode45(fun6, tspan, y0)

np.savetxt('Test16_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test16_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 17
tspan = np.linspace(10,0,25)
y0 = np.array([2])

myOptions = Options() #Create an option with default value.
myOptions.odeset('Refine',10)
myOptions.odeset('NormControl',True)
myOptions.odeset('MaxStep',1)
myOptions.odeset('InitialStep',np.array([0.1]))

sol = ode45(fun1,tspan,y0, myOptions)

np.savetxt('Test17_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test17_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 18
tspan = np.linspace(5,0,25)
y0 = np.array([0, 0.01])
A=1
B=2
arg = A,B

myOptions = Options() #Create an option with default value.
myOptions.odeset('Refine',1)
myOptions.odeset('NormControl',False)
myOptions.odeset('MaxStep',np.array([0.5]))
myOptions.odeset('InitialStep',np.array([0.1]))

sol = ode45(fun5,tspan,y0, myOptions, arg)

np.savetxt('Test18_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test18_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 19
tspan = np.linspace(3,0,25)
y0 = np.array([0, 0])

sol = ode45(fun9,tspan,y0)

np.savetxt('Test19_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test19_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 20
tspan = np.linspace(4,0,25)
y0 = np.array([0, 4, 2, 20, -np.pi/2, 2])
options = Options()
options.odeset('Mass',mass_fun8)

sol = ode45(fun8,tspan,y0,options)

np.savetxt('Test20_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test20_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 21

tspan = np.linspace(0,10,25)
y0 = np.array([2])
options = Options()
options.odeset('RelTol',np.array([1e-2]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-3)

sol = ode45(fun1, tspan, y0, options)

np.savetxt('Test21_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test21_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 22

tspan = np.linspace(0,10,25)
y0 = np.array([0])
options = Options()
options.odeset('RelTol',np.array([1e-2]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-3)

sol = ode45(fun2, tspan, y0, options)

np.savetxt('Test22_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test22_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 23

tspan = np.linspace(0,2,25)
y0 = np.array([2, 4, 6])
options = Options()
options.odeset('RelTol',np.array([1e-2]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-3)

sol = ode45(fun3, tspan, y0, options)

np.savetxt('Test23_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test23_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 24

tspan = np.linspace(0,20,25)
y0 = np.array([2, 0])
options = Options()
options.odeset('RelTol',np.array([1e-2]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-3)

sol = ode45(fun4, tspan, y0, options)

np.savetxt('Test24_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test24_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 25

A=1
B=2
arg = A,B
tspan = np.linspace(0,5,25)
y0 = np.array([0, 0.01])
options = Options()
options.odeset('RelTol',np.array([1e-2]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-3)

sol = ode45(fun5, tspan, y0, options, arg)

np.savetxt('Test25_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test25_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 26

tspan = np.linspace(0,3,25)
y0 = np.array([-3,-2,-1,0,1,2,3])
options = Options()
options.odeset('RelTol',np.array([1e-2]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-3)

sol = ode45(fun6, tspan, y0, options)

np.savetxt('Test26_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test26_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 27
tspan = np.linspace(0,10,25)
y0 = np.array([2])

options = Options()
options.odeset('RelTol',np.array([1e-2]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-3)
options.odeset('Refine',10)
options.odeset('NormControl',True)
options.odeset('MaxStep',1)
options.odeset('InitialStep',np.array([0.1]))

sol = ode45(fun1,tspan,y0, options)

np.savetxt('Test27_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test27_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 28
tspan = np.linspace(0,5,25)
y0 = np.array([0, 0.01])
A=1
B=2
arg = A,B

options = Options()
options.odeset('RelTol',np.array([1e-2]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-3)
options.odeset('Refine',1)
options.odeset('NormControl',False)
options.odeset('MaxStep',np.array([0.5]))
options.odeset('InitialStep',np.array([0.1]))

sol = ode45(fun5,tspan,y0, options, arg)

np.savetxt('Test28_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test28_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 29

tspan = np.linspace(0,3,25)
y0 = np.array([0, 0])
options = Options()
options.odeset('RelTol',np.array([1e-2]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-3)

sol = ode45(fun9,tspan,y0,options)

np.savetxt('Test29_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test29_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 30

tspan = np.linspace(0,4,25)
y0 = [0, 4, 2, 20, -np.pi/2, 2]
options = Options()
options.odeset('RelTol',np.array([1e-2]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-3)
options.odeset('Mass',mass_fun8)

sol = ode45(fun8, tspan, y0, options)

np.savetxt('Test30_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test30_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 31

tspan = np.linspace(0,10,25)
y0 = np.array([2])
options = Options()
options.odeset('RelTol',np.array([1e-6]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-8)

sol = ode45(fun1, tspan, y0, options)

np.savetxt('Test31_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test31_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 32

tspan = np.linspace(0,10,25)
y0 = np.array([0.0])
options = Options()
options.odeset('RelTol',np.array([1e-6]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-8)

sol = ode45(fun2, tspan, y0, options)

np.savetxt('Test32_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test32_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 33

tspan = np.linspace(0,2,25)
y0 = np.array([2, 4, 6])
options = Options()
options.odeset('RelTol',np.array([1e-6]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-8)

sol = ode45(fun3, tspan, y0, options)

np.savetxt('Test33_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test33_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 34

tspan = np.linspace(0,20,25)
y0 = np.array([2, 0])
options = Options()
options.odeset('RelTol',np.array([1e-6]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-8)

sol = ode45(fun4, tspan, y0, options)

np.savetxt('Test34_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test34_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 35

A=1
B=2
arg = A,B
tspan = np.linspace(0,5,25)
y0 = np.array([0, 0.01])
options = Options()
options.odeset('RelTol',np.array([1e-6]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-8)

sol = ode45(fun5, tspan, y0, options, arg)

np.savetxt('Test35_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test35_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 36

tspan = np.linspace(0,3,25)
y0 = np.array([-3,-2,-1,0,1,2,3])
options = Options()
options.odeset('RelTol',np.array([1e-6]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-8)

sol = ode45(fun6, tspan, y0, options)

np.savetxt('Test36_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test36_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 37
tspan = np.linspace(0,10,25)
y0 = np.array([2])

options = Options()
options.odeset('RelTol',np.array([1e-6]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-8)
options.odeset('Refine',10)
options.odeset('NormControl',True)
options.odeset('MaxStep',1)
options.odeset('InitialStep',np.array([0.1]))

sol = ode45(fun1,tspan,y0, options)

np.savetxt('Test37_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test37_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 38
tspan = np.linspace(0,5,25)
y0 = np.array([0, 0.01])
A=1
B=2
arg = A,B

options = Options()
options.odeset('RelTol',np.array([1e-6]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-8)
options.odeset('Refine',1)
options.odeset('NormControl',False)
options.odeset('MaxStep',np.array([0.5]))
options.odeset('InitialStep',np.array([0.1]))

sol = ode45(fun5,tspan,y0, options, arg)

np.savetxt('Test38_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test38_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 39

tspan = np.linspace(0,3,25)
y0 = np.array([0, 0])
options = Options()
options.odeset('RelTol',np.array([1e-6]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-8)

sol = ode45(fun9,tspan,y0,options)

np.savetxt('Test39_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test39_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 40

tspan = np.linspace(0,4,25)
y0 = [0, 4, 2, 20, -np.pi/2, 2]
options = Options()
options.odeset('RelTol',np.array([1e-6]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-8)
options.odeset('Mass',mass_fun8)

sol = ode45(fun8, tspan, y0, options)

np.savetxt('Test40_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test40_python_y.txt', np.transpose(sol.y), fmt='%.16e')

######################################### Test 41

tspan = np.array([0,1,2,3,4,5,6,7,8,9,10])
y0 = np.array([2])

options = Options()
options.odeset('RelTol',np.array([1e-12]))
options.odeset('AbsTol',np.ones([len(y0)])*1e-15)

sol = ode45(fun1, tspan, y0 ,options)

np.savetxt('Test41_python_t.txt', sol.t, fmt='%.16e')
np.savetxt('Test41_python_y.txt', np.transpose(sol.y), fmt='%.16e')

#########################################
#L'erreur est plus grande si le nombre de pas d'intégration sol.stats.nsteps est plus grand
#########################################