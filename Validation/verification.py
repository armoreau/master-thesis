import numpy as np
import matplotlib.pyplot as plt
from compare_data import max_rel_error

nbr_input = 35
nbr_tests = 30

Test_max_rel_error = np.zeros([nbr_input,nbr_tests])

for INPUT in range(nbr_input) :
    INPUT = INPUT+1
    Test1 = max_rel_error('Input'+ str(INPUT) + '_Test1_python.txt','Input'+ str(INPUT) + '_Test1_Matlab.txt')
    Test2 = max_rel_error('Input'+ str(INPUT) + '_Test2_python.txt','Input'+ str(INPUT) + '_Test2_Matlab.txt')
    Test3 = max_rel_error('Input'+ str(INPUT) + '_Test3_python.txt','Input'+ str(INPUT) + '_Test3_Matlab.txt')
    Test4 = max_rel_error('Input'+ str(INPUT) + '_Test4_python.txt','Input'+ str(INPUT) + '_Test4_Matlab.txt')
    Test5 = max_rel_error('Input'+ str(INPUT) + '_Test5_python.txt','Input'+ str(INPUT) + '_Test5_Matlab.txt')
    Test6 = max_rel_error('Input'+ str(INPUT) + '_Test6_python.txt','Input'+ str(INPUT) + '_Test6_Matlab.txt')
    Test7 = max_rel_error('Input'+ str(INPUT) + '_Test7_python.txt','Input'+ str(INPUT) + '_Test7_Matlab.txt')
    Test8 = max_rel_error('Input'+ str(INPUT) + '_Test8_python.txt','Input'+ str(INPUT) + '_Test8_Matlab.txt')
    Test9 = max_rel_error('Input'+ str(INPUT) + '_Test9_python.txt','Input'+ str(INPUT) + '_Test9_Matlab.txt')
    Test10 = max_rel_error('Input'+ str(INPUT) + '_Test10_python.txt','Input'+ str(INPUT) + '_Test10_Matlab.txt')
    Test11 = max_rel_error('Input'+ str(INPUT) + '_Test11_python.txt','Input'+ str(INPUT) + '_Test11_Matlab.txt')
    Test12 = max_rel_error('Input'+ str(INPUT) + '_Test12_python.txt','Input'+ str(INPUT) + '_Test12_Matlab.txt')
    Test13 = max_rel_error('Input'+ str(INPUT) + '_Test13_python.txt','Input'+ str(INPUT) + '_Test13_Matlab.txt')
    Test14 = max_rel_error('Input'+ str(INPUT) + '_Test14_python.txt','Input'+ str(INPUT) + '_Test14_Matlab.txt')
    Test15 = max_rel_error('Input'+ str(INPUT) + '_Test15_python.txt','Input'+ str(INPUT) + '_Test15_Matlab.txt')
    Test16 = max_rel_error('Input'+ str(INPUT) + '_Test16_python.txt','Input'+ str(INPUT) + '_Test16_Matlab.txt')
    Test17 = max_rel_error('Input'+ str(INPUT) + '_Test17_python.txt','Input'+ str(INPUT) + '_Test17_Matlab.txt')
    Test18 = max_rel_error('Input'+ str(INPUT) + '_Test18_python.txt','Input'+ str(INPUT) + '_Test18_Matlab.txt')
    Test19 = max_rel_error('Input'+ str(INPUT) + '_Test19_python.txt','Input'+ str(INPUT) + '_Test19_Matlab.txt')
    Test20 = max_rel_error('Input'+ str(INPUT) + '_Test20_python.txt','Input'+ str(INPUT) + '_Test20_Matlab.txt')
    Test21 = max_rel_error('Input'+ str(INPUT) + '_Test21_python.txt','Input'+ str(INPUT) + '_Test21_Matlab.txt')
    Test22 = max_rel_error('Input'+ str(INPUT) + '_Test22_python.txt','Input'+ str(INPUT) + '_Test22_Matlab.txt')
    Test23 = max_rel_error('Input'+ str(INPUT) + '_Test23_python.txt','Input'+ str(INPUT) + '_Test23_Matlab.txt')
    Test24 = max_rel_error('Input'+ str(INPUT) + '_Test24_python.txt','Input'+ str(INPUT) + '_Test24_Matlab.txt')
    Test25 = max_rel_error('Input'+ str(INPUT) + '_Test25_python.txt','Input'+ str(INPUT) + '_Test25_Matlab.txt')
    Test26 = max_rel_error('Input'+ str(INPUT) + '_Test26_python.txt','Input'+ str(INPUT) + '_Test26_Matlab.txt')
    Test27 = max_rel_error('Input'+ str(INPUT) + '_Test27_python.txt','Input'+ str(INPUT) + '_Test27_Matlab.txt')
    Test28 = max_rel_error('Input'+ str(INPUT) + '_Test28_python.txt','Input'+ str(INPUT) + '_Test28_Matlab.txt')
    Test29 = max_rel_error('Input'+ str(INPUT) + '_Test29_python.txt','Input'+ str(INPUT) + '_Test29_Matlab.txt')
    Test30 = max_rel_error('Input'+ str(INPUT) + '_Test30_python.txt','Input'+ str(INPUT) + '_Test30_Matlab.txt')

    Test_max_rel_error[INPUT-1,:] = np.array([Test1,Test2,Test3,Test4,Test5,Test6,Test7,Test8,Test9,Test10,
                       Test11,Test12,Test13,Test14,Test15,Test16,Test17,Test18,Test19,Test20,
                       Test21,Test22,Test23,Test24,Test25,Test26,Test27,Test28,Test29,Test30])

#pas grave
#Input7_test16 failure Warning: Matlab and python absh <= hmin
#Input12_Test16 idem
#input23_test16 idem
    
#Input7_test28 Data not comparable too high nstep.
#Input31_test20 Data not comparable too high nstep.

idem = 0
idem15 = 0
idem14 = 0
idem13 = 0
idem12 = 0
idem11 = 0
idem10 = 0
idem9 = 0
idem8 = 0
idem7 = 0
idem6 = 0
idem5 = 0
idem4 = 0
idem3 = 0
idem2 = 0
idem1 = 0
error =0

for i in range(nbr_input) :
    for j in range(nbr_tests) :
        if Test_max_rel_error[i,j] < 1e-50 :
            idem = idem + 1
        elif Test_max_rel_error[i,j] < 1e-15 :
            idem15 = idem15 +1
        elif Test_max_rel_error[i,j] < 1e-14 :
            idem14 = idem14 +1
        elif Test_max_rel_error[i,j] < 1e-13 :
            idem13 = idem13 +1
        elif Test_max_rel_error[i,j] < 1e-12 :
            idem12 = idem12 +1
        elif Test_max_rel_error[i,j] < 1e-11 :
            idem11 = idem11 +1
        elif Test_max_rel_error[i,j] < 1e-10 :
            idem10 = idem10 +1
        elif Test_max_rel_error[i,j] < 1e-9 :
            idem9 = idem9 +1
        elif Test_max_rel_error[i,j] < 1e-8 :
            idem8 = idem8 +1
        elif Test_max_rel_error[i,j] < 1e-7 :
            idem7 = idem7 +1
        elif Test_max_rel_error[i,j] < 1e-6 :
            idem6 = idem6 +1
        elif Test_max_rel_error[i,j] < 1e-6 :
            idem6 = idem6 +1
        elif Test_max_rel_error[i,j] < 1e-5 :
            idem5 = idem5 +1
        elif Test_max_rel_error[i,j] < 1e-4 :
            idem4 = idem4 +1
        elif Test_max_rel_error[i,j] < 1e-3 :
            idem3 = idem3 +1
        elif Test_max_rel_error[i,j] < 1e-2 :
            idem2 = idem2 +1
        elif Test_max_rel_error[i,j] < 1e-1 :
            idem1 = idem1 +1
        else :
            error = error +1
            
#Make barplot
            
#height = np.array([error, idem1, idem2, idem3, idem4, idem5, idem6, idem7, idem8, idem9, idem10, idem11, idem12, idem13, idem14, idem15, idem])
#bars = ('Erreur', '< e-1', '< e-2', '< e-3', '< e-4', '< e-5', '< e-6', '< e-7', '< e-8', '< e-9', '< e-10', '< e-11', '< e-12', '< e-13', '< e-14', '< e-15', 'Ok')

height = np.array([error, idem6, idem7, idem8, idem9, idem10, idem11, idem12, idem13, idem14, idem15, idem])
bars = ('Erreur', '< e-6', '< e-7', '< e-8', '< e-9', '< e-10', '< e-11', '< e-12', '< e-13', '< e-14', '< e-15', 'Ok')
y_pos = np.arange(len(bars))

# Create bars
plt.bar(y_pos, height)

# Create names on the x-axis
plt.xticks(y_pos, bars)

# Show graphic
plt.show()