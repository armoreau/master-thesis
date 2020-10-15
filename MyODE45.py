import numpy as np

def myODE45(fun ,t_span, y0, t_eval) :

    #Explicit Runge-kutta method order 5
    
    #Dormand-Prince coefficient
    A = np.array([
    [0, 0, 0, 0, 0, 0],
    [1/5, 0, 0, 0, 0, 0],
    [3/40, 9/40, 0, 0, 0, 0],
    [44/45, -56/15, 32/9, 0, 0, 0],
    [19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0],
    [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0],
    [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84]])

    B_1 = np.array([35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0])#next point
    B_2 = np.array([5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40])#estimation
    C = np.array([0, 1/5, 3/10, 4/5, 8/9, 1, 1])
    
    tol = 0.0001
    order = 5
    
    size_y0 = np.size(y0)
    
    output_y = np.zeros((size_y0,1))
    output_t = np.zeros((1,1))
    
    for j in range(size_y0):
        
        output_y[j,0] = y0[j]
        output_t[0] = t_span[0]
        
        t = t_span[0]
        estimation_error = 0
        
        i=0
        while t < t_span[1] :
            #step size
            if t==t_span[0] :
                h = (t_span[1]-t_span[0])/1000 #To discuss
            else :
                h = h*np.abs(tol/estimation_error)**(1/(order+1))
            if t+h > t_span[1] : #on ne depasse pas la borne d integration
                h = t_span[1] - t
            
            next_t = np.array([[t+h]])
            output_t = np.concatenate((output_t,next_t),axis=1)
            
            #7 stages
            k0 = fun(t,output_y[j,i])
            k1 = fun(t+C[1]*h,output_y[j,i]+h*(A[1,0]*k0))
            k2 = fun(t+C[2]*h,output_y[j,i]+h*(A[2,0]*k0+A[2,1]*k1))
            k3 = fun(t+C[3]*h,output_y[j,i]+h*(A[3,0]*k0+A[3,1]*k1+A[3,2]*k2))
            k4 = fun(t+C[4]*h,output_y[j,i]+h*(A[4,0]*k0+A[4,1]*k1+A[4,2]*k2+A[4,3]*k3))
            k5 = fun(t+C[5]*h,output_y[j,i]+h*(A[5,0]*k0+A[5,1]*k1+A[5,2]*k2+A[5,3]*k3+ A[5,4]*k4))
            k6 = fun(t+C[6]*h,output_y[j,i]+h*(A[6,0]*k0+A[6,1]*k1+A[6,2]*k2+A[6,3]*k3+ A[6,4]*k4+ A[6,5]*k5))

            next_y = np.array([[output_y[j,i] + h*(k0*B_1[0] + k1*B_1[1] + k2*B_1[2] + k3*B_1[3] + k4*B_1[4] + k5*B_1[5] +k6*B_1[6])]])
            output_y = np.concatenate((output_y,next_y),axis=1) #si y0 a plusieurs dim ca risque de foirer
            estimation = output_y[j,i] + h*(k0*B_2[0] + k1*B_2[1] + k2*B_2[2] + k3*B_2[3] + k4*B_2[4] + k5*B_2[5] +k6*B_2[6])
            
            estimation_error = output_y[j,i+1] - estimation
            
            i=i+1
            t=t+h
            
        if t_eval is None :
            return (output_t, output_y);
        else:
            size_t_eval = np.size(t_eval)
            new_output_y = np.zeros((size_y0,size_t_eval))
            
            k=0
            
            for j in range(size_y0):
                for i in range(size_t_eval):
                    
                    """
                    while not (output_t[j,k] <= t_eval[i] and t_eval[i] <= output_t[j,k+1]) :
                        k=k+1
                        
                    dist1 = t_eval[i]-output_t[j,k]
                    dist2 = output_t[j,k+1]-t_eval[i]
                    dist = dist1+dist2
                    
                    new_output_y[j,i] = (dist2/dist)*(output_y[j,k]+ fun(output_t[j,k],output_y[j,k])*(t_eval[i]-output_t[j,k])) + (dist1/dist)* (output_y[j,k+1] - fun(output_t[j,k+1],output_y[j,k+1])*(output_t[j,k]-t_eval[i]))"""
                    """
                    while not (t_eval[i] <= output_t[j,k+1]) :
                        k=k+1
                    
                    new_output_y[j,i] = output_y[j,k] + fun(output_t[j,k],output_y[j,k])*(t_eval[i]-output_t[j,k])"""
                                    
                    while not (output_t[j,k] <= t_eval[i] and t_eval[i] <= output_t[j,k+1]) :
                        k=k+1
                        
                    #calcul de la valeur de new_output_y[j,i] (A modifier approximation inexacte)
                    dist1 = t_eval[i]-output_t[j,k]
                    dist2 = output_t[j,k+1]-t_eval[i]
                    dist = dist1+dist2
                    
                    new_output_y[j,i] = (dist2/dist)*output_y[j,k]+ (dist1/dist)*output_y[j,k+1]

                    
            return (t_eval, new_output_y);
    
    return [output_t, output_y];

        
        
