import numpy as np
def compare_data(MatlabFileName, PythonFileName):
        
    data_Matlab_t = np.array([ float(z) for z in open(MatlabFileName, 'r').read().split() ])
    data_python_t = np.array([ float(z) for z in open(PythonFileName, 'r').read().split() ])
    
    return data_Matlab_t == data_python_t