import numpy as np
from inspect import isfunction

def myODE45(ode,tspan,y0,options = None, varargin = None) :
        
    solver_name = 'ode45'
    # Check inputs
    
    # Stats
    nsteps  = 0
    nfailed = 0
    nfevals = 0
      
    # Output
    FcnHandlesUsed  = isfunction(ode)
    if not FcnHandlesUsed :
        raise ValueError("arg1 is not an handle function")
    
    output_sol = True #(FcnHandlesUsed && (nargout==1))      # sol = odeXX(...)
    output_ty  = False #(~output_sol && (nargout > 0))  # [t,y,...] = odeXX(...)
    #There might be no output requested...
    
    # Handle solver arguments
    neq, tspan, ntspan, NEXT, t0, tfinal, tdir, y0, f0, odeArgs, odeFcn, options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, dataType = odearguments(FcnHandlesUsed, solver_name, ode, tspan, y0, options, varargin)
    nfevals = nfevals + 1
    
    #Handle the output
    
    refine = np.maximum(1,options.Refine)
    
    if ntspan > 2 :
        outputAt = 1 # output only at tspan points
    elif refine <= 1 :
        outputAt = 2 # computed points, no refinement
    else :
        outputAt = 3 # computed points, with refinement
        S = np.array(range(1,refine))/refine

    # Handle the event function
    
    # Handle the mass matrix
            
    #Non-negative solution components
    idxNonNegative = options.NonNegative
    nonNegative =  not isempty(idxNonNegative)
    thresholdNonNegative = np.abs(threshold)
    
    t = t0
    y = y0
    
    # Allocate memory if we're generating output.
    
    output_y = np.zeros((neq,1))
    output_t = np.zeros((1,1))
        
    output_y[:,0] = y0
    output_t[0] = t0
    
    #Initialize method parameters.
    POW = 1/5
    A = np.array([1/5, 3/10, 4/5, 8/9, 1, 1])
    B = np.array([
    [1/5, 3/40, 44/45, 19372/6561, 9017/3168, 35/384],
    [0, 9/40, -56/15, -25360/2187, -355/33, 0],
    [0, 0, 32/9, 64448/6561, 46732/5247, 500/1113],
    [0, 0, 0, -212/729, 49/176, 125/192],
    [0, 0, 0, 0, -5103/18656, -2187/6784],
    [0, 0, 0, 0, 0, 11/84],
    [0, 0, 0, 0, 0, 0]
    ])
    
    E = np.array([71/57600, 0, -71/16695, 71/1920, -17253/339200, 22/525, -1/40])
    f = np.zeros((neq,7))
    hmin = 16*np.finfo(float).eps
    
    if htry == None :  # Compute an initial step size h using y'(t).
        absh = np.minimum(hmax, htspan)
        if normcontrol :
            rh = (np.linalg.norm(f0) / np.maximum(normy,threshold)) / (0.8 * rtol**POW)
        else :
            rh = np.linalg.norm(f0 / np.maximum(np.abs(y0),threshold),np.inf) / (0.8 * rtol**POW)
        if absh * rh > 1 :
            absh = 1 / rh
        absh = np.maximum(absh, hmin)
    else :
        absh = np.minimum(hmax, np.maximum(hmin, htry))
    f[:,0] = f0
    #THE MAIN LOOP
    done = False
    while not done :
        
        # By default, hmin is a small number such that t+hmin is only slightly
        # different than t.  It might be 0 if t is 0.
        hmin = 16*np.finfo(float).eps
        absh = np.minimum(hmax, np.maximum(hmin, absh)) # couldn't limit absh until new hmin
        h = tdir * absh

        # Stretch the step if within 10% of tfinal-t.
        if 1.1*absh >= np.abs(tfinal - t) :
            h = tfinal - t
            absh = np.abs(h)
            done = True
   
        # LOOP FOR ADVANCING ONE STEP.
        nofailed = True  # no failed attempts
        while True :
            
            hA = h * A
            hB = h * B

            f[:,1] = feval(ode,t+hA[0],y+np.dot(f,hB[:,0]),odeArgs)
            f[:,2] = feval(ode,t+hA[1],y+np.dot(f,hB[:,1]),odeArgs)
            f[:,3] = feval(ode,t+hA[2],y+np.dot(f,hB[:,2]),odeArgs)
            f[:,4] = feval(ode,t+hA[3],y+np.dot(f,hB[:,3]),odeArgs)
            f[:,5] = feval(ode,t+hA[4],y+np.dot(f,hB[:,4]),odeArgs)
            
            tnew = t + hA[5]
            if done :
                tnew = tfinal  # Hit end point exactly.
            h = tnew - t # Purify h
            
            ynew = y + np.dot(f,hB[:,5])
            f[:,6] = feval(ode,tnew,ynew,odeArgs)
            nfevals = nfevals + 6
            
            #Estimate the error
            NNrejectStep = False
            if normcontrol :
                normynew = np.linalg.norm(ynew)
                errwt = np.maximum(np.maximum(normy,normynew),threshold)
                err = absh * (np.linalg.norm(np.dot(f,E)) / errwt)
                if nonNegative and (err <= rtol) and np.any(ynew(idxNonNegative)<0) :
                    errNN = np.linalg.norm( np.maximum(0,-ynew(idxNonNegative)) ) / errwt 
                    if errNN > rtol :
                        err = errNN
                        NNrejectStep = True    
            else :
                err = absh * (np.linalg.norm(np.dot(f,E) / np.maximum(np.maximum(np.abs(y),np.abs(ynew)),threshold),np.inf))
                if nonNegative and (err <= rtol) and np.any(ynew(idxNonNegative)<0) :
                    errNN = np.linalg.norm( np.maximum(0,-ynew(idxNonNegative)) / thresholdNonNegative, np.inf)
                    if errNN > rtol :
                        err = errNN
                        NNrejectStep = True
            
            # Accept the solution only if the weighted error is no more than the
            # tolerance rtol.  Estimate an h that will yield an error of rtol on
            # the next step or the next try at taking this step, as the case may be,
            # and use 0.8 of this value to avoid failures.
            if err > rtol :# Failed step
                nfailed = nfailed + 1            
                if absh <= hmin :
                    print("Warning:python:ode45:IntegrationTolNotMet:absh <= hmin ")
                    sol = Sol(output_t[0], output_y,nsteps,nfailed,nfevals,options)
                    return sol
          
                if nofailed :
                    nofailed = False
                    if NNrejectStep :
                        absh = np.maximum(hmin, 0.5*absh)
                    else :
                        absh = np.maximum(hmin, absh * np.maximum(0.1, 0.8*(rtol/err)**POW))
                else : 
                    absh = np.maximum(hmin, 0.5 * absh)
                h = tdir * absh
                done = False

            else : # Successful step
                NNreset_f7 = False
                if nonNegative and np.any(ynew[idxNonNegative]<0) :
                    ynew[idxNonNegative] = np.maximum(ynew[idxNonNegative],0)
                    if normcontrol :
                        normynew = np.linalg.norm(ynew)
                    NNreset_f7 = True
                break
        nsteps = nsteps + 1
             
        #GERER LES output
        if outputAt == 1 : #Evaluate only at t_span
            
            while NEXT < ntspan :
                if tspan[NEXT] == tnew :
                    to_concatanate_t = np.array([[tnew]])
                    output_t = np.concatenate((output_t,to_concatanate_t),axis=1)
                        
                    to_concatanate_y = np.transpose(np.array([ynew]))
                    output_y = np.concatenate((output_y,to_concatanate_y),axis=1)
                    
                    NEXT = NEXT+1
                elif tspan[NEXT] < tnew and t < tspan[NEXT] :
                    
                    tinterp = tspan[NEXT]
                    to_concatanate_t = np.array([[tinterp]])
                    output_t = np.concatenate((output_t,to_concatanate_t),axis=1)
                    
                    yinterp = ntrp45(tinterp,t,y,tnew,ynew,h,f)
                    to_concatanate_y = np.transpose(np.array([yinterp]))
                    output_y = np.concatenate((output_y,to_concatanate_y),axis=1)
                    
                    NEXT = NEXT+1
                elif tspan[NEXT] > tnew:
                    break
        
        else : 
            if outputAt == 3 : #Evaluate at solver steps + refined step
                for i in range(S.size) :
                    tinterp = t + h*S[i]
                    to_concatanate_t = np.array([[tinterp]])
                    output_t = np.concatenate((output_t,to_concatanate_t),axis=1)
                    
                    yinterp = ntrp45(tinterp,t,y,tnew,ynew,h,f)
                    to_concatanate_y = np.transpose(np.array([yinterp]))
                    output_y = np.concatenate((output_y,to_concatanate_y),axis=1)  
                
                to_concatanate_t = np.array([[tnew]])
                output_t = np.concatenate((output_t,to_concatanate_t),axis=1)
                
                to_concatanate_y = np.transpose(np.array([ynew]))
                output_y = np.concatenate((output_y,to_concatanate_y),axis=1)               
                
            else : #Evaluate only at solver steps
                to_concatanate_t = np.array([[tnew]])
                output_t = np.concatenate((output_t,to_concatanate_t),axis=1)
                    
                to_concatanate_y = np.transpose(np.array([ynew]))
                output_y = np.concatenate((output_y,to_concatanate_y),axis=1)
                
        
        if done :
            break
        # If there were no failures compute a new h.
        if nofailed:
        # Note that absh may shrink by 0.8, and that err may be 0.
            temp = 1.25*(err/rtol)**POW
            if temp > 0.2:
                absh = absh / temp
            else:
                absh = 5.0*absh
        # Advance the integration one step.
        t = tnew
        y = ynew
        if normcontrol:
            normy = normynew
            
        if NNreset_f7 :
        # Used f7 for unperturbed solution to interpolate.  
        # Now reset f7 to move along constraint. 
            f[:,6] = feval(ode,tnew,ynew,odeArgs)
            nfevals = nfevals + 1
        f[:,0] = f[:,6]

            
            
    sol = Sol(output_t[0], output_y,nsteps,nfailed,nfevals,options)
    return sol

def odearguments(FcnHandlesUsed, solver, ode, tspan, y0, options, extras):
    
    if FcnHandlesUsed :
        tspan = np.array(tspan)
        if tspan.size < 2 :
            raise ValueError("pyhton:odearguments:tspan.size < 2")

        htspan = abs(tspan[1] - tspan[0])
        tspan = np.array(tspan)
        ntspan = tspan.size
        t0 = tspan[0]
        NEXT = 1       # NEXT entry in tspan
        tfinal = tspan[ntspan-1]
        args = extras
        
    #else :
        #Todo

    y0 = np.array(y0)
    neq = y0.size
    
    # Test that tspan is internally consistent.
    if any(np.isnan(tspan)) :
        raise ValueError("pyhton:odearguments:TspanNaNValues")
    if t0 == tfinal :
        raise ValueError("pyhton:odearguments:TspanEndpointsNotDistinct")
    
    if tfinal > t0 :
        tdir = 1
    else :
        tdir = -1
    #tdir = np.sign(tfinal - t0)
    #print(tdir)
        
    if any( tdir*np.diff(tspan) <= 0 ) :
        raise ValueError("pyhton:odearguments:TspanNotMonotonic")
    
    f0 = feval(ode,t0,y0,args) #Regler le prob de args
    
    if options == None :     
        options = Options(neq,tfinal-t0) #Use default values
    #if options.AbsTol == None :
        options.AbsTol = 1e-6*np.ones(neq)
    if options.MaxStep == None:
        options.MaxStep = np.abs(0.1*(tfinal-t0))
        
    rtol = np.array(options.RelTol)
    if (rtol.size != 1 or rtol <= 0) :
        raise ValueError("pyhton:odearguments:RelTolNotPosScalar")
    if rtol < 100*np.finfo(float).eps :
        rtol = 100*np.finfo(float).eps
        
    atol = np.array(options.AbsTol)
    #if any(atol <= 0) :
        #raise ValueError("python:odearguments:AbsTolNotPos")
        
    normcontrol = options.NormControl
    if normcontrol :
        if atol.size != 1 :
            raise ValueError("python:odearguments:NonScalarAbsTol")
        normy = np.linalg.norm(y0)
    else :
        if ((atol.size != 1) and (atol.size != neq)) :
            raise ValueError("python:odearguments:SizeAbsTol")
        atol = np.array(atol)
        normy = None
            
    threshold = atol/rtol
        
    hmax = options.MaxStep
    if hmax <= 0 :
        raise ValueError("python:odearguments:MaxStepLEzero")
        
    htry = options.InitialStep
    if htry != None :
        if (not isempty(htry) and (htry <= 0)) :
            raise ValueError("python:odearguments:InitialStepLEzero")
        
    odeFcn = ode
    dataType = None
    
    return neq, tspan, ntspan, NEXT, t0, tfinal, tdir, y0, f0, args, odeFcn, options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, dataType
        

def ntrp45(tinterp,t,y,tnew,ynew,h,f): #use to interpolate the solution at t_span point
    BI = np.array([
    [1, -183/64, 37/12, -145/128],
    [0, 0 ,0 ,0],
    [0, 1500/371, -1000/159, 1000/371],
    [0, -125/32, 125/12, -375/64],
    [0, 9477/3392, -729/106, 25515/6784],
    [0, -11/7, 11/3, -55/28],
    [0, 3/2, -4, 5/2]])
 
    s = ((tinterp - t) / h)*np.ones(4)
    
    yinterp = y + np.dot(np.dot(f,(h*BI)),np.cumprod(s))
    return yinterp

def feval(fun,t,y,args) :
    if args == None :
        return np.array(fun(t,y))
    else :
        return np.array(fun(t,y,*args))

def isempty(x):
    if x.size == 0:
        return True
    return False

class Options :
        
    def __init__(self,neq=None,length_tspan=None):
        self.RelTol = np.array([1e-3])
        if neq == None :
            self.AbsTol = None
        else :
            self.AbsTol = 1e-6*np.ones(neq)
        self.NormControl = False
        if length_tspan == None :
            self.MaxStep = None
        else :
            self.MaxStep = np.abs(0.1*(length_tspan))
        self.InitialStep = None
        self.Refine = 4
        self.NonNegative = np.array([])
        
    def odeset(self,string,value) :
        if string == 'RelTol' :
            self.RelTol = value
        elif string == 'AbsTol':
            self.AbsTol = value
        elif string == 'NormControl' :
            self.NormControl = value
        elif string == 'MaxStep' :
            self.MaxStep = value
        elif string == 'InitialStep' :
            self.InitialStep = value
        elif string == 'Refine' :
            self.Refine = value
        elif string == 'NonNegative' :
            self.NonNegative = value
        else :
            print("Warning:odeset:OptionsNotFoundIgnoreg")
        
              
        
class Sol :
    def __init__(self,t,y,nsteps,nfailed,nfevals,options):
        self.t = t
        self.y = y
        self.nsteps = nsteps
        self.nfailed = nfailed
        self.nfevals = nfevals
        self.options = options
