import numpy as np

class Options :
        
    def __init__(self,neq=None,length_tspan=None):
        self.RelTol = np.array([1e-3])
        if neq == None :
            self.AbsTol = None
        else :
            self.AbsTol = 1e-6*np.ones(neq)
        self.NormControl = False
        if length_tspan == None :
            self.MaxStep = None
        else :
            self.MaxStep = np.abs(0.1*(length_tspan))
        self.InitialStep = None
        self.Refine = 4
        self.NonNegative = np.array([])
        
    def odeset(self,string,value) :
        if string == 'RelTol' :
            self.RelTol = value
        elif string == 'AbsTol':
            self.AbsTol = value
        elif string == 'NormControl' :
            self.NormControl = value
        elif string == 'MaxStep' :
            self.MaxStep = value
        elif string == 'InitialStep' :
            self.InitialStep = value
        elif string == 'Refine' :
            self.Refine = value
        elif string == 'NonNegative' :
            self.NonNegative = value
        else :
            print("Warning:odeset:OptionsNotFoundIgnoreg")