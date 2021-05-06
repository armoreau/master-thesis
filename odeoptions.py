import numpy as np

class Odeoptions :
        
    def __init__(self,neq=None,length_tspan=None):
        
        if neq is None :
            self.AbsTol = None
        else :
            self.AbsTol = 1e-6*np.ones(neq)
        if length_tspan is None :
            self.MaxStep = None
        else :
            self.MaxStep = np.abs(0.1*(length_tspan))
            
        self.NormControl = False  
        self.RelTol = np.array([1e-3])
        self.InitialStep = None
        self.Refine = 4
        self.NonNegative = np.array([])
        self.Events = None
        self.Mass = None
        self.MStateDependence = 'weak'
        
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
        elif string == 'Events' :
            self.Events = value
        elif string == 'Mass' :
            self.Mass = value
        elif string == 'MStateDependence' :
            self.MStateDependence = value
        else :
            raise Exception("python:odearguments:odeset:OptionsNotFound")