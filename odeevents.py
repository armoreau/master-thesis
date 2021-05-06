import numpy as np
from feval import feval

def odeevents(FcnHandlesUsed,ode,t0,y0,options,extras):
    
    haveeventfun = 0   # no Events function
    eventArgs = None
    eventValue = None
    teout = np.array([])
    yeout = np.array([])
    ieout = np.array([])
    
    eventFcn = options.Events
    if eventFcn is None :
      return haveeventfun, eventFcn, eventArgs, eventValue, teout, yeout, ieout
    
    if FcnHandlesUsed :    # function handles used 
      haveeventfun = 1   # there is an Events function
      eventArgs = extras
      [eventValue, trash1, trash2] = feval(eventFcn,t0,y0,eventArgs)
    
    return haveeventfun, eventFcn, eventArgs ,eventValue, teout, yeout, ieout