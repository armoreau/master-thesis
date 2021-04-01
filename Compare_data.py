import numpy as np
import math

#arrondi un vecteur avec digit chiffre significatif
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