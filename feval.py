import numpy as np

def feval(fun,t,y,args) :
    if args is None :
        return np.array(fun(t,y))

    else :
        return np.array(fun(t,y,*args))