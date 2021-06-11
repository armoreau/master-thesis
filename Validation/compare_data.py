import numpy as np
import math

#rounded vector with digit significant digit
def signif(x, digit):
    toRet = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] == 0:
            toRet[i] = 0
        else :
            toRet[i] = round(x[i], digit - int(math.floor(math.log10(abs(x[i])))) - 1)
    return toRet

def compare_data(MatlabFileName, PythonFileName, digit):
        
    data_Matlab = np.array([ float(z) for z in open(MatlabFileName, 'r').read().split() ])
    data_python = np.array([ float(z) for z in open(PythonFileName, 'r').read().split() ])
    
    rounded_data_Matlab = signif(data_Matlab, digit)
    rounded_data_python = signif(data_python, digit)
    
    Boolean_vector = rounded_data_Matlab == rounded_data_python
    return Boolean_vector.all()

def max_rel_error(MatlabFileName, PythonFileName):
    
    data_Matlab = np.array([ float(z) for z in open(MatlabFileName, 'r').read().split() ])
    data_python = np.array([ float(z) for z in open(PythonFileName, 'r').read().split() ])
    
    if(len(data_Matlab) != len(data_python)) :
        return 1000 #Key for not comparable data
    
    abs_diff = np.abs(data_Matlab - data_python)
    rel_diff = abs_diff
    
    length = len(abs_diff)
    for i in range(length) :
        if (np.abs(data_Matlab[i]) > 1) :
            rel_diff[i] = abs_diff[i]/np.abs(data_Matlab[i])
            
    if len(rel_diff) == 0:
        return 0
    else :
        return np.max(rel_diff)