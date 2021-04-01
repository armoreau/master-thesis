import numpy as np
import matplotlib.pyplot as plt
from Options import Options
from ode45 import ode45

def events(t,y):
    value = np.array([y[0]])
    isterminal = np.array([1])
    direction = np.array([-1])
    return [value,isterminal,direction]

def dydt(t,y):
    return np.array([y[1],-9.8])

tstart = 0
tfinal = 30
y0 = np.array([0.0, 20.0])
options = Options()
options.odeset('Refine',100) #More beautiful
options.odeset('Events',events)

tout = np.array([tstart])
yout = np.array([[y0[0]],[y0[1]]])

teout = np.array([0])
yeout = np.array([[0],[0]])
ieout = np.array([0])

for i in range(10) :
    
    tspan=np.array([tstart,tfinal]) 
    res = ode45(dydt,tspan,y0,options)
    nt = len(res.t)
    
    teout=np.concatenate((teout,res.te))
    yeout=np.concatenate((yeout,res.ye),axis=1)
    ieout=np.concatenate((ieout,res.ie))
    
    tout=np.concatenate((tout,res.t[1:len(res.t)]))
    yout=np.concatenate((yout,res.y[:,1:len(res.t)]),axis=1)

    y0[0] = 0
    y0[1] = -0.9*res.y[1,nt-1]
    
    options.odeset('InitialStep',np.array([res.t[nt-1]-res.t[nt-4-1]]))
    options.odeset('MaxStep',res.t[nt-1]-res.t[0])
    
    tstart = res.t[nt-1]
    
#PLOT
fig = plt.figure()
plt.title('Ball trajectory and the events')
plt.xlabel('time')
plt.ylabel('height')
plt.plot(tout,yout[0],label='position')
plt.plot(teout, yeout[0], 'ro')
plt.show()
