import numpy as np

from feval import feval
from ntrp45 import ntrp45
from isempty import isempty


def odezero(ntrpfun, eventfun, eventargs, v, t, y, tnew, ynew, t0, h, f, idxNonNegative) :
    
    # Initialize.
    tol = 128*np.maximum(np.finfo(float).eps,np.finfo(float).eps)
    tol = np.minimum(tol, np.abs(tnew - t))
    tout = np.array([])
    yout = np.array([[],[]])
    iout = np.array([])
    
    tdir = np.sign(tnew - t)
    stop = 0
    rmin = np.finfo(float).eps

    # Set up tL, tR, yL, yR, vL, vR, isterminal and direction.
    tL = t
    yL = y
    vL = v
    [vnew,isterminal,direction] = feval(eventfun,tnew,ynew,eventargs)
    if isempty(direction) :
        direction = np.zeros(len(vnew))   # zeros crossings in any direction
    tR = tnew
    yR = ynew
    vR = vnew

    # Initialize ttry so that we won't extrapolate if vL or vR is zero.
    ttry = tR

    # Find all events before tnew or the first terminal event.
    while True :
        lastmoved = 0
        while True :
            # Events of interest shouldn't have disappeared, but new ones might
            # be found in other elements of the v vector.
            indzc = [i for i in range(len(direction)) if np.sign(vR[i])!=np.sign(vL[i]) and direction[i]*(vR[i]-vL[i])>=0]
            if isempty(indzc) :
                if lastmoved != 0 :
                    raise Exception('ode45:odezero:LostEvent')
                else :
                    return tout,yout,iout,vnew,stop
            
            # Check if the time interval is too short to continue looking.
            delta = tR - tL
            if np.abs(delta) <= tol:
                break
            
            if (tL == t) and any([vL[index]==0 and vR[index]!=0 for index in indzc]) :
                ttry = tL + tdir*0.5*tol
             
            else :
                #Compute Regula Falsi change, using leftmost possibility.
                change = 1
                for j in indzc :
                    # If vL or vR is zero, try using old ttry to extrapolate.
                    if vL[j]== 0:
                        if (tdir*ttry > tdir*tR) and (vtry[j] != vR[j]) :
                            maybe = 1.0 - vR[j] * (ttry-tR) / ((vtry[j]-vR[j]) * delta)
                            if (maybe < 0) or (maybe > 1) :
                                maybe = 0.5
                        else:
                            maybe = 0.5
                    elif vR[j] == 0.0:
                        if (tdir*ttry < tdir*tL) and (vtry[j] != vL[j]) :
                            maybe = vL[j] * (tL-ttry) / ((vtry[j]-vL) * delta)
                            if (maybe < 0) or (maybe > 1):
                                maybe = 0.5
                        else:
                            maybe = 0.5
                    else:
                        maybe = -vL[j] / (vR[j] - vL[j]) #Note vR(j) != vL(j).
                    if maybe < change:
                        change = maybe
                change = change * np.abs(delta)
                
                # Enforce minimum and maximum change.
                change = np.maximum(0.5*tol, np.minimum(change, np.abs(delta) - 0.5*tol))
                ttry = tL + tdir * change
                    
            # Compute vtry.
            ytry, trash1 = ntrp45(ttry,t,y,tnew,ynew,h,f,idxNonNegative)
            [vtry, trash2, trash3] = feval(eventfun,ttry,ytry,eventargs)
              
            # Check for any crossings between tL and ttry.
            indzc=[i for i in range(len(direction)) if np.sign(vtry[i])!=np.sign(vL[i]) and direction[i]*(vtry[i]-vL[i])>=0]
            if (not isempty(indzc)):
                # Move right end of bracket leftward, remembering the old value.
                tswap = tR
                tR = ttry
                ttry = tswap
                yswap = yR
                yR = ytry
                ytry = yswap
                vswap = vR
                vR = vtry
                vtry = vswap
                # Illinois method.  If we've moved leftward twice, halve
                # vL so we'll move closer next time.
                if lastmoved == 2:
                    # Watch out for underflow and signs disappearing.
                    maybe = 0.5 * vL
                    i = [j for j in range(len(maybe)) if np.abs(maybe[j] >= rmin)]
                    for temp in i:
                        vL[temp] = maybe[temp]
                        
                lastmoved = 2
            else:
                #Move left end of bracket rightward, remembering the old value.
                tswap = tL
                tL = ttry
                ttry = tswap
                yswap = yL
                yL = ytry
                ytry = yswap
                vswap = vL
                vL = vtry
                vtry = vswap
                # Illinois method.  If we've moved rightward twice, halve
                # vR so we'll move closer next time.
                if lastmoved == 1:
                    # Watch out for underflow and signs disappearing.
                    maybe = 0.5 * vR
                    i = [j for j in range(len(maybe)) if np.abs(maybe[j] >= rmin)]
                    for temp in i:
                        vR[temp] = maybe[temp]
                        
                lastmoved = 1

        j = np.ones([len(indzc)])
        add_tout=np.array([tR for i in j])
        add_yout=np.tile(np.transpose(np.array([yR[:,0]])),len(indzc))
        add_iout=np.transpose(np.array([indzc]))
        if isempty(tout) :
            tout=add_tout
            yout=add_yout
            iout=add_iout
        else :
            tout=np.concatenate((tout,add_tout))
            yout=np.concatenate((yout,add_yout),axis=1)
            iout=np.concatenate((iout,add_iout))
                
        if any([isterminal[i] for i in indzc]):
            if tL != t0:
                stop = 1
            break
        elif np.abs(tnew - tR) <= tol :
            #  We're not going to find events closer than tol.
            break
        else :
            # Shift bracket rightward from [tL tR] to [tR+0.5*tol tnew].
            ttry = tR
            ytry = yR
            vtry = vR
            tL = tR + tdir*0.5*tol          
            yL, trash1 = ntrp45(tL,t,y,tnew,ynew,h,f,idxNonNegative)
            [vL, trash2, trash3] = feval(eventfun,tL,yL,eventargs)            
            tR = tnew
            yR = ynew
            vR = vnew
             
    return tout,yout,iout,vnew,stop
