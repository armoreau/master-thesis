import numpy as np

def feval(fun,t,y,args) :
    if args is None :
        if y is None :
            return np.array(fun(t))
        else :
            return np.array(fun(t,y))
    else :
        if y is None :
            return np.array(fun(t,*args))
        else :
            return np.array(fun(t,y,*args))
        
