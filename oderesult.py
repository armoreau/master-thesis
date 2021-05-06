class Oderesult :
    def __init__(self,solver,extdata,t,y,stats,teout,yeout,ieout):
        self.solver = solver
        self.extdata = extdata
        self.t = t
        self.y = y
        self.stats = stats
        self.te = teout
        self.ye = yeout
        self.ie = ieout

class Extdata :
    def __init__(self,odefun,options,varargin):
        self.odefun = odefun
        self.options = options
        self.varargin = varargin

class Stats :
    def __init__(self,nsteps,nfailed,nfevals):
        self.nsteps = nsteps
        self.nfailed = nfailed
        self.nfevals = nfevals
        