import numpy as np
from odeoptions import Odeoptions
from feval import feval

def odearguments(FcnHandlesUsed, solver, ode, tspan, y0, options, extras):
    
    if FcnHandlesUsed :
        tspan = np.array(tspan)
        if tspan.size < 2 :
            raise Exception("pyhton:odearguments:tspan.size < 2")

        htspan = np.abs(tspan[1] - tspan[0])
        tspan = np.array(tspan)
        ntspan = tspan.size
        t0 = tspan[0]
        NEXT = 1       # NEXT entry in tspan
        tfinal = tspan[ntspan-1]
        args = extras

    y0 = np.array(y0)
    neq = len(y0)
    
    # Test that tspan is internally consistent.
    if any(np.isnan(tspan)) :
        raise Exception("pyhton:odearguments:TspanNaNValues")
    if t0 == tfinal :
        raise Exception("pyhton:odearguments:TspanEndpointsNotDistinct")
    
    if tfinal > t0 :
        tdir = 1
    else :
        tdir = -1
        
    if any( tdir*np.diff(tspan) <= 0 ) :
        raise Exception("pyhton:odearguments:TspanNotMonotonic")
    
    f0 = feval(ode,t0,y0,args)
    
    if options is None :     
        options = Odeoptions() #Use default values
        
    if options.MaxStep is None :
        options.MaxStep = np.abs(0.1*(tfinal-t0))
        
    rtol = np.array([options.RelTol])
    if (len(rtol) != 1 or rtol <= 0) :
        raise Exception("pyhton:odearguments:RelTolNotPosScalar")
    if rtol < 100*np.finfo(float).eps :
        rtol = 100*np.finfo(float).eps
    
    atol = options.AbsTol
    if isinstance(atol,list):
        atol = np.array(atol)
    else :
        atol = np.array([atol])
        
    if any(atol <= 0) :
        raise Exception("python:odearguments:AbsTolNotPos")
        
    normcontrol = options.NormControl
    if normcontrol :
        if len(atol) != 1 :
            raise Exception("python:odearguments:NonScalarAbsTol")
        normy = np.linalg.norm(y0)
    else :
        if ((len(atol) != 1) and (len(atol) != neq)) :
            raise Exception("python:odearguments:SizeAbsTol")
        normy = None
            
    threshold = atol/rtol
        
    hmax = np.array([options.MaxStep])
    if hmax <= 0 :
        raise Exception("python:odearguments:MaxStepLEzero")
        
    htry = options.InitialStep
    if htry is not None :
        if htry <= 0 :
            raise Exception("python:odearguments:InitialStepLEzero")
        
    odeFcn = ode
    dataType = 'float64'
    
    return neq, tspan, ntspan, NEXT, t0, tfinal, tdir, y0, f0, args, odeFcn, options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, dataType