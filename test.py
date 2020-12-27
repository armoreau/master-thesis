import numpy as np
from inspect import isfunction
import math
import MyODE45


#def myfun(ode,tspan,y0,options = None, varargin = None) :
#    
#    args = varargin
#    print(args)
#    
#    return 45
#
#def eq3(t, y):
#    return -np.cos(1/y[0])/(y**2)
#
#class Option:
#
#    def __init__(self):
#        self.rTol = 10e-4
#
#option = Option()

m = np.array([5])

c= np.array([2])
a= np.array([3])
b= np.array([1])

arg = []

string = 'coucou'
print(string == 'cOucou')