import numpy as np   
        
class Sol :
    def __init__(self,t,y,nsteps,nfailed,nfevals,options,teout, yeout, ieout):
        self.t = t
        self.y = y
        self.nsteps = nsteps
        self.nfailed = nfailed
        self.nfevals = nfevals
        self.options = options
        self.teout = teout
        self.yeout = yeout
        self.ieout = ieout