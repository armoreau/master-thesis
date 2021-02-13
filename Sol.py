import numpy as np   
        
class Sol :
    def __init__(self,t,y,nsteps,nfailed,nfevals,options):
        self.t = t
        self.y = y
        self.nsteps = nsteps
        self.nfailed = nfailed
        self.nfevals = nfevals
        self.options = options