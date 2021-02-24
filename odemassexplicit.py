import numpy as np
from feval import feval

def ExplicitSolverHandleMass(t,y,odeFcn,massFcn,varargin) :
    A = feval(massFcn,t,y,varargin)
    b = feval(odeFcn,t,y,varargin)
    yp = np.linalg.lstsq(A, b,rcond=None)[0]

    return yp

def odemassexplicit(FcnHandlesUsed,massType,odeFcn,odeArgs,massFcn,massM) :
    if FcnHandlesUsed :
        new_odeFcn = ExplicitSolverHandleMass
        odeArgs = (odeFcn,massFcn,odeArgs)
        return new_odeFcn,odeArgs
    else : #ode-file use TODO
        return odeFcn,odeArgs