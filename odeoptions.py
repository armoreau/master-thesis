class Odeoptions :
        
    def __init__(self,length_tspan=None):
         
        self.RelTol = 1e-3
        self.AbsTol = 1e-6
        self.NormControl = False 
        self.MaxStep = None 
        self.InitialStep = None
        self.Refine = 4
        self.NonNegative = []
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