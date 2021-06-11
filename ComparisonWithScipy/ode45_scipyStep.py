import numpy as np
from inspect import isfunction

#Add parent folder to the path. Code taken from https://codeolives.com/2020/01/10/python-reference-module-in-parent-directory/
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from oderesult import Oderesult,Stats,Extdata
from feval import feval
from isempty import isempty
from ntrp45 import ntrp45
from odearguments import odearguments
from odeevents import odeevents
from odezero import odezero
from odemass import odemass
from odemassexplicit import odemassexplicit
from odenonnegative import odenonnegative

def ode45_scipyStep(ode,tspan,y0,scipyStep,options = None, varargin = None) :
        
    solver_name = 'ode45'
    # Check inputs
    
    # Stats
    nsteps  = 0
    nfailed = 0
    nfevals = 0
      
    # Output
    FcnHandlesUsed  = isfunction(ode)
    if not FcnHandlesUsed :
        raise ValueError("ode is not an handle function")
    
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
    haveEventFcn,eventFcn,eventArgs,valt,teout,yeout,ieout = odeevents(FcnHandlesUsed,odeFcn,t0,y0,options,varargin)
    
    # Handle the mass matrix
    Mtype, M, Mfun = odemass(FcnHandlesUsed,odeFcn,t0,y0,options,varargin)
    if Mtype > 0 :
        #check if matrix is singular and raise an arror
        
        # Incorporate the mass matrix into odeFcn and odeArgs.
        odeFcn,odeArgs = odemassexplicit(FcnHandlesUsed,Mtype,odeFcn,odeArgs,Mfun,M)
        f0 = feval(odeFcn,t0,y0,odeArgs)
        nfevals = nfevals + 1
            
    #Non-negative solution components
    idxNonNegative = options.NonNegative
    nonNegative =  not isempty(idxNonNegative)
    #thresholdNonNegative = np.abs(threshold)
    if nonNegative :
        odeFcn,thresholdNonNegative,odeArgs = odenonnegative(odeFcn,y0,threshold,idxNonNegative,odeArgs)
        #odeArgs = (argSup,odeArgs)
        f0 = feval(odeFcn,t0,y0,odeArgs)
        nfevals = nfevals + 1
    
    t = t0
    y = y0
    
    # Allocate memory if we're generating output.
    
    if outputAt == 1 :
        output_y = np.zeros((neq,ntspan),dtype=dataType)
        output_t = np.zeros(ntspan,dtype=dataType)
    else :
        output_y = np.zeros((neq,1),dtype=dataType)
        output_t = np.zeros(1,dtype=dataType)
        
    output_y[:,0] = y0
    output_t[0] = t0
    
    #Initialize method parameters.
    POW = 1/5
    A = np.array([1/5, 3/10, 4/5, 8/9, 1, 1],dtype=dataType)
    B = np.array([
    [1/5, 3/40, 44/45, 19372/6561, 9017/3168, 35/384],
    [0, 9/40, -56/15, -25360/2187, -355/33, 0],
    [0, 0, 32/9, 64448/6561, 46732/5247, 500/1113],
    [0, 0, 0, -212/729, 49/176, 125/192],
    [0, 0, 0, 0, -5103/18656, -2187/6784],
    [0, 0, 0, 0, 0, 11/84],
    [0, 0, 0, 0, 0, 0]
    ],dtype=dataType)
    
    E = np.array([71/57600, 0, -71/16695, 71/1920, -17253/339200, 22/525, -1/40],dtype=dataType)
    f = np.zeros((neq,7),dtype=dataType)
    hmin = 16*np.finfo(float(t)).eps
    
    if htry is None :  # Compute an initial step size h using y'(t).
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
    stop = False
    
    plop = 0
    while not done :
        
        h = scipyStep[plop]
        plop = plop +1
        absh = h
        # Stretch the step if within 10% of tfinal-t.
        if plop >= len(scipyStep) :
            done = True
   
        # LOOP FOR ADVANCING ONE STEP.
        nofailed = True  # no failed attempts
        
        while True :
                     
            hA = h * A
            hB = h * B

            f[:,1] = feval(odeFcn,t+hA[0],y+np.dot(f,hB[:,0]),odeArgs)
            f[:,2] = feval(odeFcn,t+hA[1],y+np.dot(f,hB[:,1]),odeArgs)
            f[:,3] = feval(odeFcn,t+hA[2],y+np.dot(f,hB[:,2]),odeArgs)
            f[:,4] = feval(odeFcn,t+hA[3],y+np.dot(f,hB[:,3]),odeArgs)
            f[:,5] = feval(odeFcn,t+hA[4],y+np.dot(f,hB[:,4]),odeArgs)
            
            tnew = t + hA[5]
            if done :
                tnew = tfinal  # Hit end point exactly.
            h = tnew - t # Purify h
            
            ynew = y + np.dot(f,hB[:,5])
            f[:,6] = feval(odeFcn,tnew,ynew,odeArgs)
            nfevals = nfevals + 6

            #else : # Successful step
            NNreset_f7 = False
            if nonNegative and np.any(ynew[idxNonNegative]<0) :
                ynew[idxNonNegative] = np.maximum(ynew[idxNonNegative],0)
                if normcontrol :
                    normynew = np.linalg.norm(ynew)
                NNreset_f7 = True
            break
        nsteps = nsteps + 1
             
        if haveEventFcn :
            te,ye,ie,valt,stop=odezero(None,eventFcn,eventArgs,valt,t,y,tnew,ynew,t0,h,f,idxNonNegative)
            
            if not isempty(te) :
                
                teout=np.append(teout,te)
                if isempty(yeout) :
                    yeout=ye
                else:
                    yeout=np.append(yeout,ye,axis=1)
                ieout=np.append(ieout,ie)
                
                if stop :
                    taux = t + (te[-1] - t)*A
                    trash, f[:,1:7]=ntrp45(taux,t,y,None,None,h,f,idxNonNegative)
                    tnew = te[-1]
                    ynew = ye[:,-1]
                    h = tnew - t
                    done = True

        #GERER LES output
        if outputAt == 1 : #Evaluate only at t_span
            
            if tdir == 1 : #tspan is increasing
                while NEXT < ntspan :
                    if tspan[NEXT] == tnew :
                        
                        output_t[NEXT] = tnew     
                        output_y[:,NEXT] = ynew                       
                        NEXT = NEXT+1
                        
                    elif tspan[NEXT] < tnew and t < tspan[NEXT] :
                        
                        first_indice = NEXT
                        NEXT = NEXT+1
                        while tspan[NEXT] < tnew :
                            NEXT = NEXT+1
                        last_indice = NEXT
                        
                        tinterp = tspan[first_indice:last_indice]
                        output_t[first_indice:last_indice] = tinterp
                        
                        yinterp = ntrp45(tinterp,t,y,tnew,ynew,h,f,idxNonNegative,Need_ypinterp = False)
                        output_y[:,first_indice:last_indice] = yinterp
                        
                    elif stop == True :
                        output_t = np.append(output_t,[tnew],axis=0)
                        
                        to_concatenate_y = np.transpose(np.array([ynew]))
                        output_y = np.append(output_y,to_concatenate_y,axis=1)
                        
                        break
                        
                    elif tspan[NEXT] > tnew:
                        break
                    
            elif tdir == -1 : #tspan is decreasing
                while NEXT < ntspan :
                    if tspan[NEXT] == tnew :
                        
                        output_t[NEXT] = tnew
                        output_y[:,NEXT] = ynew
                        NEXT = NEXT+1
                        
                    elif tspan[NEXT] > tnew and t > tspan[NEXT] :
                        
                        first_indice = NEXT
                        NEXT = NEXT+1
                        while tspan[NEXT] > tnew :
                            NEXT = NEXT+1
                        last_indice = NEXT
                        
                        tinterp = tspan[first_indice:last_indice]
                        output_t[first_indice:last_indice] = tinterp
                        
                        yinterp = ntrp45(tinterp,t,y,tnew,ynew,h,f,idxNonNegative,Need_ypinterp = False)
                        output_y[:,first_indice:last_indice] = yinterp
                        
                    elif stop == True : # DEBUG
                        output_t = np.append(output_t,[tnew],axis=0)
                        
                        to_concatenate_y = np.transpose(np.array([ynew]))
                        output_y = np.append(output_y,to_concatenate_y,axis=1)
                        
                        break
                    
                    elif tspan[NEXT] < tnew:
                        break
        
        else : 
            if outputAt == 3 : #Evaluate at solver steps + refined step
                    
                tinterp = t + h*S               
                output_t = np.append(output_t,tinterp,axis=0)
                
                yinterp = ntrp45(tinterp,t,y,tnew,ynew,h,f,idxNonNegative,Need_ypinterp = False)
                output_y = np.append(output_y,yinterp,axis=1)
                
                output_t = np.append(output_t,[tnew],axis=0)
                
                to_concatenate_y = np.transpose(np.array([ynew]))
                
                output_y = np.append(output_y,to_concatenate_y,axis=1)               
                
            else : #Evaluate only at solver steps
                output_t = np.append(output_t,[tnew],axis=0)
                    
                to_concatenate_y = np.transpose(np.array([ynew]))
                output_y = np.append(output_y,to_concatenate_y,axis=1)
                
        
        if done :
            break
        
        # Advance the integration one step.
        t = tnew
        y = ynew
        if normcontrol:
            normy = normynew
            
        if NNreset_f7 :
        # Used f7 for unperturbed solution to interpolate.  
        # Now reset f7 to move along constraint. 
            f[:,6] = feval(odeFcn,tnew,ynew,odeArgs)
            nfevals = nfevals + 1
        f[:,0] = f[:,6]
        
    extdata = Extdata(odeFcn,options,odeArgs)
    stats = Stats(nsteps,nfailed,nfevals)
    oderesult = Oderesult(solver_name,extdata,output_t,output_y,stats,teout,yeout,ieout)
    return oderesult