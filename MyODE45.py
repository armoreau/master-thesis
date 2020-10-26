import numpy as np

# Multiply steps computed from asymptotic behaviour of errors by this.
SAFETY = 0.9

MIN_FACTOR = 0.2  # Minimum allowed decrease in a step size.
MAX_FACTOR = 10  # Maximum allowed increase in a step size.

ERROR_EXPONENT = -1 / (4 + 1) # 4 est l'ordre de l'erreur d'estimation

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

    B = np.array([35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0])#next point
    E = np.array([5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40])#estimation
    C = np.array([0, 1/5, 3/10, 4/5, 8/9, 1, 1])
    
    tol = 0.0001
    order = 4
    
    size_y0 = np.size(y0)
    
    output_y = np.zeros((size_y0,1))
    output_t = np.zeros((1,1))
        
    output_y[:,0] = y0
    output_t[0] = t_span[0]
        
    t = t_span[0]
    y=y0
        
    h=0
    while t < t_span[1] :
            #step size
        if t==t_span[0] :
            direction = 1
            f0 = np.array(fun(t,y0))
            h = choose_first_step(fun, t, y0, f0, direction, order)
            
        max_step = t_span[1] - t_span[0]
        rtol = 1e-3#100 * EPS
        atol = 1e-6#4 * EPS

        min_step = 10 * np.abs(np.nextafter(t, direction * np.inf) - t)
            
        if h > max_step:
            h = max_step
        elif h < min_step:
            h = min_step

        step_accepted = False
        step_rejected = False
            
        while not step_accepted:

            h = h*direction
            next_t = t + h

            if direction * (next_t - t_span[1]) > 0:
                next_t = t_span[1]

            h = next_t - t
            h = np.abs(h)

            next_y, estimation_error = step(fun,t,y,h,A,B,C,E)
            scale = atol + np.maximum(np.abs(y), np.abs(next_y)) * rtol
                           
            error_norm = norm(estimation_error / scale)

            if error_norm < 1:
                if error_norm == 0:
                    factor = MAX_FACTOR
                else:
                    factor = min(MAX_FACTOR, SAFETY * error_norm ** ERROR_EXPONENT)

                if step_rejected:
                    factor = min(1, factor)

                h *= factor

                step_accepted = True
            else:
                h *= max(MIN_FACTOR, SAFETY * error_norm ** ERROR_EXPONENT)
                step_rejected = True
            
        #next_y = step(fun,t,y,h,A,B,C)
        #max_estimation_error = step_error(fun,t,y,h,A,B,C,E)
        
        to_concatanate_y = np.transpose(np.array([next_y]))
        output_y = np.concatenate((output_y,to_concatanate_y),axis=1)
        
        #next_t = t+h
        to_concatanate_t = np.array([[next_t]])
        output_t = np.concatenate((output_t,to_concatanate_t),axis=1)
            
        y = next_y
        t = next_t
            
    return (output_t[0], output_y);  
        
def choose_first_step(fun, t0, y0, f0, direction, order):
    EPS = np.finfo(float).eps
    rtol = 1e-3#100 * EPS
    atol = 1e-6#4 * EPS

    scale = atol + np.abs(y0) * rtol
    d0 = norm(y0 / scale)
    d1 = norm(f0 / scale)
    if d0 < 1e-5 or d1 < 1e-5:
        h0 = 1e-6
    else:
        h0 = 0.01 * d0 / d1

    y1 = y0 + h0 * direction * f0
    f1 = fun(t0 + h0 * direction, y1)
    d2 = norm((f1 - f0) / scale) / h0

    if d1 <= 1e-15 and d2 <= 1e-15:
        h1 = max(1e-6, h0 * 1e-3)
    else:
        h1 = (0.01 / max(d1, d2)) ** (1 / (order + 1))

    return min(100 * h0, h1)
    
def norm(x):
    """Compute RMS norm."""
    return np.linalg.norm(x) / x.size ** 0.5

def step(fun,t,y,h,A,B,C,E): #a ameliore (pas besoin de calculer 2 fois les coeff K)
    
    k0 = np.array(fun(t,y))
    k1 = np.array(fun(t+C[1]*h,y+h*(A[1,0]*k0)))
    k2 = np.array(fun(t+C[2]*h,y+h*(A[2,0]*k0+A[2,1]*k1)))
    k3 = np.array(fun(t+C[3]*h,y+h*(A[3,0]*k0+A[3,1]*k1+A[3,2]*k2)))
    k4 = np.array(fun(t+C[4]*h,y+h*(A[4,0]*k0+A[4,1]*k1+A[4,2]*k2+A[4,3]*k3)))
    k5 = np.array(fun(t+C[5]*h,y+h*(A[5,0]*k0+A[5,1]*k1+A[5,2]*k2+A[5,3]*k3+ A[5,4]*k4)))
    k6 = np.array(fun(t+C[6]*h,y+h*(A[6,0]*k0+A[6,1]*k1+A[6,2]*k2+A[6,3]*k3+ A[6,4]*k4+ A[6,5]*k5)))
    
    next_y = y + h*(k0*B[0] + k1*B[1] + k2*B[2] + k3*B[3] + k4*B[4] + k5*B[5] +k6*B[6])
    estimation = y + h*(k0*E[0] + k1*E[1] + k2*E[2] + k3*E[3] + k4*E[4] + k5*E[5] +k6*E[6])
    estimation_error = next_y - estimation

    
    return next_y, estimation_error
