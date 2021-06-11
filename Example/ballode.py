import numpy as np
import matplotlib.pyplot as plt

#Add parent folder to the path. Code taken from https://codeolives.com/2020/01/10/python-reference-module-in-parent-directory/
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from odeoptions import Odeoptions
from ode45 import ode45

def events(t,y):
    value = [y[0]]
    isterminal = [1]
    direction = [-1]
    return [value,isterminal,direction]

def dydt(t,y):
    return [y[1],-9.8]

tstart = 0
tfinal = 30
y0 = [0.0, 20.0]
options = Odeoptions()
options.odeset('Refine',10)
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
    
    options.odeset('InitialStep',res.t[nt-1]-res.t[nt-4-1])
    options.odeset('MaxStep',res.t[nt-1]-res.t[0])
    
    tstart = res.t[nt-1]
    
#PLOT
fig = plt.figure()
plt.title('Ball trajectory and the events')
plt.xlabel('time [s]')
plt.ylabel('height [m]')
plt.plot(tout,yout[0],label='position')
plt.plot(teout, yeout[0], 'ro')
plt.show()