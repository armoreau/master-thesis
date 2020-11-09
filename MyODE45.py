import numpy as np

def myODE45(fun ,t_span, y0) :

    nfevals = 0
    nfailed = 0
    nsteps = 0
    
    size_y0 = np.size(y0)
    
    output_y = np.zeros((size_y0,1))
    output_t = np.zeros((1,1))
        
    output_y[:,0] = y0
    output_t[0] = t_span[0]
        
    t = t_span[0]
    y = y0
    f0 = np.array(fun(t,y0))
    
    #Compute variable like matlab :
    hmin,htspan,hmax,atol,rtol,threshold,normcontrol,POW,tdir = odearguments(t_span,size_y0)
    tfinal = t_span[t_span.size-1] #a mettre dans odearguments
    idxNonNegative = np.array([])
    nonNegative =  not isempty(idxNonNegative)
    thresholdNonNegative = np.abs(threshold)
    
    if t_span.size > 2 :
        RequestedPoints = True
        NEXT=1
    else :
        RequestedPoints = False
        refine = 4
        if refine > 1 :
            RefinedSteps = True #Evaluate sol at refined steps
            S = np.array(range(1,refine))/refine
        else :
            RefinedSteps = False #Evaluate sol at solver steps
    neq = size_y0

    POW = 1/5

    #Explicit Runge-kutta method order 5
    
    #Dormand-Prince coefficient
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

    #B = np.array([35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0])#next point
    #E = np.array([5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40])#estimation
    

    f = np.zeros((neq,7))
    hmin = 16*np.finfo(float).eps
    absh = choose_first_step_matlab(hmin, hmax, htspan, y0, f0, threshold, rtol, normcontrol, POW)
    f[:,0] = f0
    
    #main loop
    done = False
    while not done :
        
        hmin = 16*np.finfo(float).eps
        absh = np.minimum(hmax, np.maximum(hmin, absh))
        h = tdir * absh

        if 1.1*absh >= np.abs(tfinal - t) :
            h = tfinal - t
            absh = np.abs(h)
            done = True
   
        # LOOP FOR ADVANCING ONE STEP.
        nofailed = True  
        while True :
            
            #ynew, estimation_error = step(fun,t,y,h,A,B,C,E)
            hA = h * A
            hB = h * B
            
            f[:,1] = np.array(fun(t+hA[0],y+np.dot(f,hB[:,0])))
            f[:,2] = np.array(fun(t+hA[1],y+np.dot(f,hB[:,1])))
            f[:,3] = np.array(fun(t+hA[2],y+np.dot(f,hB[:,2])))
            f[:,4] = np.array(fun(t+hA[3],y+np.dot(f,hB[:,3])))
            f[:,5] = np.array(fun(t+hA[4],y+np.dot(f,hB[:,4])))
            #f[:,0] = np.array(fun(t+C[6]*h,y+h*(A[6,0]*k0+A[6,1]*k1+A[6,2]*k2+A[6,3]*k3+ A[6,4]*k4+ A[6,5]*k5)))
            
            tnew = t + hA[5]
            if done :
                tnew = tfinal  
            h = tnew - t
            
            ynew = y + np.dot(f,hB[:,5])
            f[:,6] = fun(tnew,ynew)
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
                    print("Error absh <= hmin")
                    return
          
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

            else :# Successful step
                NNreset_f7 = False
                if nonNegative and np.any(ynew[idxNonNegative]<0) :
                    ynew[idxNonNegative] = np.maximum(ynew[idxNonNegative],0)
                    if normcontrol :
                        normynew = norm(ynew)
                    NNreset_f7 = True                 
                break
        nsteps = nsteps + 1
             
        #GERER LES output
        if RequestedPoints : #Evaluate only at t_span
            
            while NEXT < t_span.size :
                if t_span[NEXT] == tnew :
                    to_concatanate_t = np.array([[tnew]])
                    output_t = np.concatenate((output_t,to_concatanate_t),axis=1)
                        
                    to_concatanate_y = np.transpose(np.array([ynew]))
                    output_y = np.concatenate((output_y,to_concatanate_y),axis=1)
                    
                    NEXT = NEXT+1
                elif t_span[NEXT] < tnew and t < t_span[NEXT] :
                    
                    tinterp = t_span[NEXT]
                    to_concatanate_t = np.array([[tinterp]])
                    output_t = np.concatenate((output_t,to_concatanate_t),axis=1)
                    
                    yinterp = ntrp45(tinterp,t,y,tnew,ynew,h,f)
                    to_concatanate_y = np.transpose(np.array([yinterp]))
                    output_y = np.concatenate((output_y,to_concatanate_y),axis=1)
                    
                    NEXT = NEXT+1
                elif t_span[NEXT] > tnew:
                    break
        
        else :
            if RefinedSteps : #Evaluate at solver steps + refined step
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
            f[:,6] = fun(tnew,ynew)
            nfevals = nfevals + 1
        f[:,0] = f[:,6]
            
    return (output_t[0], output_y)

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


def choose_first_step_matlab(hmin, hmax, htspan, y0, f0, threshold, rtol, normcontrol, POW):
    absh = np.minimum(hmax, htspan)
    if normcontrol :
        rh = (np.linalg.norm(f0) / np.maximum(np.linalg.norm(y0),threshold)) / (0.8 * rtol**POW)
    else :
        rh = np.linalg.norm(f0 / np.maximum(np.abs(y0),threshold),np.inf) / (0.8 * rtol**POW)
    if absh * rh > 1 :
        absh = 1 / rh
    absh = np.maximum(absh, hmin)
    return absh

def odearguments(t_span,size_y0):
    hmin = np.finfo(float).eps
    htspan = t_span[t_span.size-1]-t_span[0]
    hmax = 0.1*(htspan)
    atol = 1e-6*np.ones(size_y0)
    rtol = 1e-3
    if rtol < 100*np.finfo(float).eps :
        rtol = 100*np.finfo(float).eps
    threshold = atol/rtol
    normcontrol = False
    POW = 1/5
    tdir = 1
    return hmin,htspan,hmax,atol,rtol,threshold,normcontrol,POW,tdir
def rms_norm(x):
    """Compute RMS norm."""
    return np.linalg.norm(x) / x.size ** 0.5

def norm(x):
    """Compute norm."""
    return np.linalg.norm(x)

def isempty(x):
    if x.size == 0:
        return True
    return False
