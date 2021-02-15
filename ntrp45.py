import numpy as np

def ntrp45(tinterp,t,y,tnew,ynew,h,f,idxNonNegative): #use to interpolate the solution at t_span point
    BI = np.array([
    [1, -183/64, 37/12, -145/128],
    [0, 0 ,0 ,0],
    [0, 1500/371, -1000/159, 1000/371],
    [0, -125/32, 125/12, -375/64],
    [0, 9477/3392, -729/106, 25515/6784],
    [0, -11/7, 11/3, -55/28],
    [0, 3/2, -4, 5/2]])
 
    #s = ((tinterp - t) / h)*np.ones(4)
    #yinterp = y + np.dot(np.dot(f,(h*BI)),np.cumprod(s))
    #return yinterp
    
    if type(tinterp) is not type(np.array([])):
        tinterp=np.array([tinterp])
    
    s = ((tinterp - t) / h)
    yinterp = np.transpose(np.tile(y,(len(tinterp),1))) + np.dot(np.dot(f,(h*BI)),np.cumprod(np.tile(s,(4,1)),axis=0))
    
    ncumprod = np.concatenate(([np.ones(len(s))],np.cumprod(np.array([2*s,(3/2)*s,(4/3)*s]),axis=0)),axis=0)
    ypinterp = np.dot(np.dot(f,BI),ncumprod)
    
    return yinterp, ypinterp