import numpy as np
import random

def Random_tspan() :
    first = random.uniform(-5, 10)
    last = random.uniform(first, 15)
    
    if random.random() < 0.5 : #50%
        res = np.array([first, last])
    else :
        npoint = random.randint(5, 20)
        if random.random() < 0.5 : #equidistant
            res = np.linspace(first,last,npoint)
        else :
            res = np.zeros(npoint)
            for i in range(npoint) :
                if i == 0 :
                    res[i] = first + random.random()
                else :
                    res[i] = res[i-1] + random.random()
    if random.random() < 0.25 : #decresing order
        temp = np.zeros(len(res))
        for i in range(len(res)) :
            temp[i] = res[len(res)-1-i]
        return temp
    else :
        return res

def Random_y0(n) :
    res = np.zeros(n)
    a = -8
    b = 8
    for i in range(n) :
        res[i] = random.uniform(a, b)
    return res

def Random_varargin() :
    res = np.zeros(3)
    for i in range(3) :
        res[i] = random.uniform(-5, 5)
    return res

####Helper for option

def Random_RelTol():
    if random.random() < 0.5 : #50%
        return 0 #default value use
    return 10**(-random.uniform(2, 10))
def Random_AbsTol():
    if random.random() < 0.5 : #50%
        return 0 #default value use
    return 10**(-random.uniform(2, 10))

def Random_Maxstep():
    if random.random() < 0.5 : #50%
        return 0 #default value use
    return 10**(random.uniform(-2, 2))

def Random_InitialStep():
    if random.random() < 0.5 : #50%
        return 0 #default value use
    return 10**(random.uniform(-6, 1))

def Random_Normcontrol():
    if random.random() < 0.5 : #50%
        return 0 #default value use
    if random.random() < 0.5 : #50%
        return True
    else :
        return False

def Random_Refine() :
    if random.random() < 0.5 : #50%
        return 0 #default value use
    return  random.randint(1, 50)

def Random_NonNegative() :
    if random.random() < 0.5 : #50%
        return 0 #default value use
    res = np.zeros(3)
    for i in range(3) :
        if random.random() < 0.5 : #50%
            res[i] = True
        else :
            res[i] = False
    return res

def Random_options(): 
    if random.random() < 0.5 : #50%
        return 0 #default value use
    else :
        return Random_RelTol(),Random_AbsTol(),Random_Maxstep(),Random_InitialStep(),Random_Normcontrol(),Random_Refine(),Random_NonNegative()

class Random_input :
    def __init__(self):
        self.tspan = Random_tspan()
        self.y01 = Random_y0(1)
        self.y02 = Random_y0(2)
        self.y03 = Random_y0(3)
        self.varargin = Random_varargin()
        self.options = Random_options()