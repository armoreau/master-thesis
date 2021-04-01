import numpy as np
from inspect import isfunction
from Sol import Sol,Stats,Extdata
from feval import feval
from isempty import isempty
from ntrp45 import ntrp45
from odearguments import odearguments
from odeevents import odeevents
from odezero import odezero
from odemass import odemass
from odemassexplicit import odemassexplicit
from odenonnegative import odenonnegative

def ode45(ode,tspan,y0,options = None, varargin = None) :
        
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
    haveEventFcn,eventFcn,eventArgs,valt,teout,yeout,ieout = odeevents(FcnHandlesUsed,odeFcn,t0,y0,options,varargin)
    
    # Handle the mass matrix
    Mtype, M, Mfun, massArgs, dMoptions = odemass(FcnHandlesUsed,odeFcn,t0,y0,options,varargin)
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
    while not done :
        # By default, hmin is a small number such that t+hmin is only slightly
        # different than t.  It might be 0 if t is 0.
        hmin = 16*np.finfo(float(t)).eps
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
            
            #Estimate the error
            NNrejectStep = False
            if normcontrol :
                normynew = np.linalg.norm(ynew)
                errwt = np.maximum(np.maximum(normy,normynew),threshold)
                err = absh * (np.linalg.norm(np.dot(f,E)) / errwt)
                if nonNegative and (err <= rtol) and np.any(ynew[idxNonNegative]<0) :
                    errNN = np.linalg.norm( np.maximum(0,-ynew[idxNonNegative]) ) / errwt 
                    if errNN > rtol :
                        err = errNN
                        NNrejectStep = True    
            else :
                err = absh * (np.linalg.norm(np.dot(f,E) / np.maximum(np.maximum(np.abs(y),np.abs(ynew)),threshold),np.inf))
                if nonNegative and (err <= rtol) and np.any(ynew[idxNonNegative]<0) :
                    errNN = np.linalg.norm( np.maximum(0,-ynew[idxNonNegative]) / thresholdNonNegative, np.inf)
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
                    extdata = Extdata(odeFcn,options,odeArgs)
                    stats = Stats(nsteps,nfailed,nfevals)
                    sol = Sol(solver_name,extdata,output_t,output_y,stats,teout,yeout,ieout)
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
            f[:,6] = feval(odeFcn,tnew,ynew,odeArgs)
            nfevals = nfevals + 1
        f[:,0] = f[:,6]
        
    extdata = Extdata(odeFcn,options,odeArgs)
    stats = Stats(nsteps,nfailed,nfevals)
    sol = Sol(solver_name,extdata,output_t,output_y,stats,teout,yeout,ieout)
    return sol



