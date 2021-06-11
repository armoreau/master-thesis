import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.integrate import solve_ivp

#Add parent folder to the path. Code taken from https://codeolives.com/2020/01/10/python-reference-module-in-parent-directory/
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from ode45 import ode45
from odeoptions import Odeoptions

######################################### Function test
### First order function
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

#### second order or more function

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

#### additionnal argument function
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

#### event function

def f26(t,y):
    dydt = np.zeros(2)
    dydt[0] = y[1]
    dydt[1] = -9.8
    return dydt

def events_f26(t,y):
    value = [y[0]]
    isterminal = [1]
    direction = [-1]
    return [value,isterminal,direction]

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
    direction = [1, -1]
    return [value,isterminal,direction]

def f28(t,y) :
    dydt = np.zeros(3)
    dydt[0] = y[1]
    dydt[1] = 1
    dydt[2] = 2*y[1]
    return dydt

def mass_f28(t,y) :
    M = np.zeros([3,3])
    M[0,0] = 1
    M[1,1] = y[1]+1
    M[1,2] = 0
    M[2,2] = 1
    return M

def f29(t,y) :
    dydt = np.zeros(2)
    dydt[0] = t
    dydt[1] = 3*t+2
    return dydt

def mass_f29(t,y) :
    M = np.zeros([2,2])
    M[0,0] = 1
    M[0,1] = 0
    M[1,0] = t
    M[1,1] = 1
    return M

def f30(t,y) :
    dydt = 2*y +t
    return dydt

def mass_f30(t,y) :
    M = np.zeros([1,1])
    M[0,0] = 5
    return M

######################################### Test

def compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT) :
    
    arg = None
    ### len(y0)=1
    y0 = y0_1
    sol = ode45(f1, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '1' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f2, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '2' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f3, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '3' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f4, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '4' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f5, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '5' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f6, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '6' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f7, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '7' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f8, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '8' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f9, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '9' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f10, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '10' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f11, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '11' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f12, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '12' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f13, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '13' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f14, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '14' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f15, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '15' + '_python.txt', np.transpose(sol.y), fmt='%.16e')

    ### len(y0)=2
    y0 = y0_2
    sol = ode45(f16, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '16' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f17, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '17' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f18, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '18' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f19, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '19' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f20, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '20' + '_python.txt', np.transpose(sol.y), fmt='%.16e')

    ### len(y0)=3
    y0 = y0_3
    sol = ode45(f21, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '21' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f22, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '22' + '_python.txt', np.transpose(sol.y), fmt='%.16e')

    ### arg function
    arg = A,B
    y0 = y0_2
    sol = ode45(f23, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '23' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    
    arg = A,B,C
    y0 = y0_1
    sol = ode45(f24, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '24' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    sol = ode45(f25, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '25' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    
    arg = None
    
    ### event function
    y0 = [0.0, 20.0]
    opts.odeset('Events',events_f26)
    sol = ode45(f26, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '26' + '_python.txt', np.transpose(sol.ye), fmt='%.16e')
    
    y0 = [1.2, 0, 0, -1.04935750983031990726]
    opts.odeset('Events',events_f27)
    sol = ode45(f27, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '27' + '_python.txt', np.transpose(sol.ye), fmt='%.16e')
    
    opts.odeset('Events',None)
    
    ### mass function
    y0 = y0_3
    opts.odeset('Mass',mass_f28)
    opts.odeset('MStateDependence','weak')
    sol = ode45(f28, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '28' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    
    y0 = y0_2
    opts.odeset('Mass',mass_f29)
    opts.odeset('MStateDependence','none')
    sol = ode45(f29, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '29' + '_python.txt', np.transpose(sol.y), fmt='%.16e')

    y0 = y0_1
    opts.odeset('Mass',mass_f30)
    opts.odeset('MStateDependence','weak')
    sol = ode45(f30, tspan, y0, opts, arg)
    np.savetxt('Input' + str(INPUT) + '_Test' + '30' + '_python.txt', np.transpose(sol.y), fmt='%.16e')
    
    opts.odeset('Mass',None)
    
######################################### perform test with random input
#### INPUT 1
INPUT = 1
tspan = [ 4.60789919,  5.06415666,  5.87811384,  6.12675463,  6.16305106,
        6.85479862,  7.12163822,  8.0260298 ,  8.86827706,  9.22210906,
       10.13889398, 10.63595004, 10.94778138, 10.99829212, 11.1745613 ,
       11.93196793, 12.90115849, 13.68223449, 14.1079588 , 14.85334451]
y0_1 = [4.53897499]
y0_2 = [1.20051259, 1.17666892]
y0_3 = [ 0.45052893,  7.03086018, -1.50526611]
A = -2.83009697
B = -4.97673659
C = -3.50342477
opts = Odeoptions()
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

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
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 3
INPUT = INPUT + 1
tspan =  [-1.41580028, 5.64584431, 7.95851231]
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
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

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
opts.odeset('Refine',4)
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

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
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 6
INPUT = INPUT + 1
tspan = [10.20618698, 10.8141241 , 11.53241461, 11.69010707, 11.86686873,
       12.57758618, 12.99199835, 13.94074236, 14.25437019]
y0_1 = [-2.51324011]
y0_2 = [-6.18446334,  7.34621481]
y0_3 = [-4.52846858,  2.46373529,  2.13908458]
A = 0.97802556
B = -2.07916483
C = -0.90394038
opts = Odeoptions()
opts.odeset('InitialStep',0.006940920311305178)
opts.odeset('NormControl',False)
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 7
INPUT = INPUT + 1
tspan = [14.34676291, -0.48332446]
y0_1 = [1.45812793]
y0_2 = [-0.23465422, -3.00879572]
y0_3 = [ 5.20014711,  3.49909892, -2.43879495]
A = -2.52625172
B = -3.07288758
C = 2.71940342
opts = Odeoptions()
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 8
INPUT = INPUT + 1
opts = Odeoptions()
tspan = [ 0.5835929 , 11.41241183]
y0_1 = [-6.16053854]
y0_2 = [4.62962776, 6.26526089]
y0_3 = [-6.78523373, -7.61812467, -6.62845498]
A = 4.82456837
B = 4.87407354
C = -0.51263804
opts = Odeoptions()
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

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
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 10
INPUT = INPUT + 1
opts = Odeoptions()
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
opts.odeset('Refine',22)
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 11
INPUT = INPUT + 1
tspan = [6.06630453, 6.85595327, 7.00007251, 7.41078496, 7.85695321]
y0_1 = [-6.64751259]
y0_2 = [-3.5624895 , -1.92005363]
y0_3 = [ 1.19939346,  5.62223293, -6.47515068]
A = -3.6500774
B = 1.32572362
C = 0.16942695
opts = Odeoptions()
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 12
INPUT = INPUT + 1
tspan = [10.33600232,  9.83453822]
y0_1 = [-4.2694042]
y0_2 = [7.47046585, 4.87288227]
y0_3 = [-6.23101346, -7.49398607, -2.95299713]
A = -0.16891248
B = 3.4227567
C = 0.33005933
opts = Odeoptions()
opts.odeset('AbsTol',0.00028328386769520237)
opts.odeset('MaxStep',0.41009555370457595)
opts.odeset('InitialStep',0.03300773645455483)
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 13
INPUT = INPUT + 1 
tspan = [3.73942976,  4.55831047,  5.37719119,  6.19607191,  7.01495262,
        7.83383334,  8.65271406,  9.47159477, 10.29047549, 11.1093562 ,
       11.92823692]
y0_1 = [0.0485922]
y0_2 = [-2.63842309,  3.56543058]
y0_3 = [-7.42857628,  7.52946529, -1.44846586]
A = 4.60762159
B = 2.90300831
C = 1.71290846
opts = Odeoptions()
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 14
INPUT = INPUT + 1 
tspan = [0.23664909, 0.40719876, 0.57774842, 0.74829809, 0.91884776,
       1.08939742, 1.25994709, 1.43049675, 1.60104642, 1.77159608,
       1.94214575, 2.11269542, 2.28324508, 2.45379475, 2.62434441,
       2.79489408, 2.96544375, 3.13599341, 3.30654308, 3.47709274]
y0_1 = [5.46724578]
y0_2 = [-7.31855781,  2.52616015]
y0_3 = [ 4.50310439, -6.63021287,  7.36660554]
A = 2.17565639
B = 2.29407597
C = -4.95602666
opts = Odeoptions()
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

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
opts.odeset('NonNegative',[0])
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 16
INPUT = INPUT + 1 
tspan = [ 3.56324789,  4.14352114,  4.7237944 ,  5.30406765,  5.88434091,
        6.46461416,  7.04488741,  7.62516067,  8.20543392,  8.78570717,
        9.36598043,  9.94625368, 10.52652694, 11.10680019, 11.68707344,
       12.2673467 ]
y0_1 = [-2.26297594]
y0_2 = [ 2.02273623, -6.45190008]
y0_3 = [-0.10547361,  1.9468299 , -7.04713625]
A = 2.36580799
B = -3.06154501
C = -2.35234947
opts = Odeoptions()
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

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
opts.odeset('MaxStep',5.244320432179101)
opts.odeset('NormControl',False)
opts.odeset('Refine',6)

compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 18
INPUT = INPUT + 1 
tspan = [2.0505969 , 2.15719485, 2.2637928 , 2.37039075, 2.4769887 ,
       2.58358665, 2.6901846 , 2.79678255, 2.9033805 , 3.00997844,
       3.11657639, 3.22317434, 3.32977229, 3.43637024, 3.54296819,
       3.64956614, 3.75616409, 3.86276204, 3.96935999, 4.07595794]
y0_1 = [7.12694417]
y0_2 = [ 3.56975967, -2.07817288]
y0_3 = [ 7.98461333, -0.59689978, -4.61550035]
A = -3.47937704
B = -3.0336046
C = 2.73797931
opts = Odeoptions()
opts.odeset('RelTol',0.007654990276229565)
opts.odeset('InitialStep',0.0008513865364431812)
opts.odeset('Refine',22)
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

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
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 20
INPUT = INPUT + 1 
tspan = [ 9.46232509,  9.73067385,  9.99902261, 10.26737138, 10.53572014,
       10.8040689 , 11.07241766, 11.34076642, 11.60911518, 11.87746394,
       12.1458127 ]
y0_1 = [5.7059878]
y0_2 = [-3.07163682,  5.98669745]
y0_3 = [ 2.54796126, -7.00795209,  1.20623878]
A = 3.44085258
B = -4.57285375
C = -3.82489802
opts = Odeoptions()
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 21
INPUT = INPUT + 1 
tspan = [3.35966743, 4.03978865, 4.71990986, 5.40003108, 6.0801523 ,
       6.76027352]
y0_1 = [7.73557472]
y0_2 = [ 2.96211931, 0.53208961]
y0_3 = [5.89125929, 6.10059695,  3.19421586]
A = 3.81728708
B = 2.10322666
C = 3.58740611
opts = Odeoptions()
opts.odeset('AbsTol',0.00016540018580299407)
opts.odeset('MaxStep',0.6327093376488327)
opts.odeset('InitialStep',0.08911214687593723)
opts.odeset('NormControl',False)
opts.odeset('NonNegative',[0])
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

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
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

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
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 24
INPUT = INPUT + 1 
tspan = [6.89384038, 6.47275395, 6.05166752, 5.63058109, 5.20949466,
       4.78840823, 4.3673218 , 3.94623537, 3.52514894]
y0_1 = [0.06673227]
y0_2 = [ 4.06547772, -5.77667385]
y0_3 = [-6.09485915,  2.83317466,  2.84113184]
A = -2.76688306
B = 2.19180466
C = -4.79997074
opts = Odeoptions()
opts.odeset('Refine',1)
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 25
INPUT = INPUT + 1 
tspan = [ 7.57582624,  7.61049506,  8.54204564,  8.79538296,  9.706019  ,
       10.41823361, 10.76703301, 10.78236865, 10.85191597, 11.23227003,
       11.93507158, 12.09226992, 12.12520625, 12.82869073, 13.77347972,
       13.89512411, 14.67831353, 15.31323283, 16.15936357]
y0_1 = [1.49109626]
y0_2 = [-4.56886996, -4.18827044]
y0_3 = [ 7.00094593, -1.76127202, -1.00658487]
A = 1.978158
B = -3.53996799
C = 1.24866477
opts = Odeoptions()
opts.odeset('Refine',3)
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 26
INPUT = INPUT + 1 
tspan = [ 9.29247421, 11.20154662]
y0_1 = [6.42056626]
y0_2 = [ 3.02751907, 7.23128958]
y0_3 = [6.8328569 ,  7.95585027,  5.25184783]
A = 3.6811953
B = 1.09621923
C = 1.53552888
opts = Odeoptions()
opts.odeset('MaxStep',0.2701987578019222)
opts.odeset('InitialStep',0.0018545051990431387)
opts.odeset('NonNegative',[0])
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

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
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

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
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 29
INPUT = INPUT + 1 
tspan = [4.40592519, 4.74395191, 5.08197864, 5.42000536, 5.75803209,
       6.09605881, 6.43408554, 6.77211226]
y0_1 = [-1.88405234]
y0_2 = [3.41725858, 2.64320194]
y0_3 = [ 3.32851184,  2.69142253, -5.39137005]
A = -4.89867833
B = 4.9015022
C = -1.20383235
opts = Odeoptions()
opts.odeset('RelTol',1.4648154972812207e-08)
opts.odeset('Refine',24)
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 30
INPUT = INPUT + 1 
tspan = [-2.29121844, -0.87546871,  0.54028101,  1.95603074,  3.37178046,
        4.78753019,  6.20327991,  7.61902964,  9.03477936]
y0_1 = [1.7769748]
y0_2 = [-7.99034351,  1.84091706]
y0_3 = [3.97775301, 6.53477137, 6.98213957]
A = 3.48665123
B = -1.78313751
C = -1.70433862
opts = Odeoptions()
opts.odeset('AbsTol',9.284105257511226e-2)
opts.odeset('MaxStep',19.31479184306004)
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

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
opts.odeset('MaxStep',0.0335244417456466)
opts.odeset('InitialStep',1.062281191880226e-05)
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 32
INPUT = INPUT + 1 
tspan = [0.41345667, 1.22961299, 1.71983443, 2.17209572, 2.24521564]
y0_1 = [-2.8119527]
y0_2 = [2.83055902, 0.41888443]
y0_3 = [7.74890026, -0.01720299,  7.74763048]
A = 4.21513872
B = 0.61825742
C = -2.4283297
opts = Odeoptions()
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

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
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 34
INPUT = INPUT + 1 
tspan = [ 7.25734918,  8.02703545,  8.79672172,  9.56640799, 10.33609426,
       11.10578053, 11.8754668 ]
y0_1 = [-1.4258659]
y0_2 = [-5.23934392, -1.76197689]
y0_3 = [-5.65019124, -5.02152462,  7.71309572]
A = -2.10258408
B = -1.58306468
C = 0.67760786
opts = Odeoptions()
opts.odeset('RelTol',0.0007149302009860351)
opts.odeset('MaxStep',7.5247957216589585)
opts.odeset('NormControl',False)
opts.odeset('Refine',4)
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)

#### INPUT 35
INPUT = INPUT + 1 
tspan = [ 9.69755071, 10.11348351, 10.52941631, 10.94534911, 11.36128191,
       11.77721471, 12.19314751, 12.60908031, 13.02501311, 13.44094591,
       13.85687871, 14.27281151]
y0_1 = [4.91590426]
y0_2 = [-5.53783951, -1.2282395 ]
y0_3 = [-5.35875016, -1.18568651,  0.73626608]
A = -4.68940883
B = 4.87638203
C = 2.92662521
opts = Odeoptions()
compute_tests(tspan,y0_1,y0_2,y0_3,opts,A,B,C,INPUT)