import numpy as np
from feval import feval

def local_odeFcn_nonnegative(t,y,ode,idxNonNegative,varargin = None) :
    yp = feval(ode,t,y,varargin)
    ndx = [i for i in range(len(idxNonNegative)) if y[i] <= 0]
    for i in ndx :
        yp[i] = np.maximum(yp[i],0)
    return yp


def odenonnegative(ode,y0,threshold,idxNonNegative,odeArgs) :
    neq = len(y0)
    thresholdNonNegative = np.array([])
    if any( (idxNonNegative < 0) or (idxNonNegative > neq) ) :
        raise Exception("python:odenonnegative:NonNegativeIndicesInvalid")

    if any(y0[idxNonNegative] < 0) :
        raise Exception("python:odenonnegative:NonNegativeViolatedAtT0")
  
    if len(threshold) == 1 :
        thresholdNonNegative = threshold[np.zeros(len(idxNonNegative),dtype=int)]
    else :
        thresholdNonNegative = threshold[idxNonNegative]

    odeFcn = local_odeFcn_nonnegative
    odeArgs = (ode,idxNonNegative,odeArgs)
    return odeFcn,thresholdNonNegative,odeArgs