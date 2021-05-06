import numpy as np
from feval import feval

def odemass(FcnHandlesUsed,ode,t0,y0,options,extras):
    
    massType = 0 
    massFcn = None
    massArgs = None
    massM = np.eye(len(y0))  
    dMoptions = None    # options for odenumjac computing d(M(t,y)*v)/dy
     
    Moption = options.Mass
    if Moption is None :
        return massType, massM, massFcn #, massArgs, dMoptions
    elif not callable(Moption):
        massType = 1
        massM = Moption
        return massType, massM, massFcn
    else : #Moption is a matrix function
        massFcn = Moption
        massArgs = extras
        Mstdep = options.MStateDependence
        if Mstdep == 'none': # time-dependent only
            massType = 2
            massM = feval(massFcn,t0,None,massArgs)
        elif Mstdep == 'weak': # state-dependent
            massType = 3
            massM = feval(massFcn,t0,y0,massArgs)
        else:
            raise Exception("python:odemass:MStateDependenceMassType")
            
        return massType, massM, massFcn #, massArgs, dMoptions
