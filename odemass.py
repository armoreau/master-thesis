import numpy as np
from feval import feval
from Options import Options

def odemass(FcnHandlesUsed,ode,t0,y0,options,extras):
    
    massType = 0 
    massFcn = None
    massArgs = None
    massM = np.eye(len(y0))  
    dMoptions = None    # options for odenumjac computing d(M(t,y)*v)/dy
    
    if FcnHandlesUsed  :   # function handles used    
        Moption = options.Mass
        if Moption is None :
            return massType, massM, massFcn, massArgs, dMoptions
        else : #Moption is a matrix function
            massType = 1 #different from matlab
            massFcn = Moption
            massArgs = extras
            massM = feval(massFcn,t0,y0,massArgs)
            
            return massType, massM, massFcn, massArgs, dMoptions
    else : #ode-file use TODO
        return massType, massM, massFcn, massArgs, dMoptions
