import numpy as np

def rkf45(ode_function, tspan, y0, tolerance=1.e-8):
    # Fehlberg coefficients
    a = [0, 1/4, 3/8, 12/13, 1, 1/2]
    b = [
        [0, 0, 0, 0, 0],
        [1/4, 0, 0, 0, 0],
        [3/32, 9/32, 0, 0, 0],
        [1932/2197, -7200/2197, 7296/2197, 0, 0],
        [439/216, -8, 3680/513, -845/4104, 0],
        [-8/27, 2, -3544/2565, 1859/4104, -11/40]
    ]
    c4 = [25/216, 0, 1408/2565, 2197/4104, -1/5, 0]
    c5 = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55]
    
    t0, tf = tspan
    t = t0
    y = np.array(y0, dtype=np.float64)
    tout = [t]
    yout = [y.copy()]
    
    h = (tf - t0) / 100  # Assumed initial time step

    while t < tf:
        hmin = 16 * np.finfo(float).eps
        ti = t
        yi = y.copy()
        f = np.zeros((len(y0), 6))
        
        # Evaluate the time derivative(s) at six points within the current interval
        for i in range(6):
            t_inner = ti + a[i] * h
            y_inner = yi.copy()
            for j in range(i):
                y_inner += h * b[i][j] * f[:, j]
            f[:, i] = ode_function(t_inner, y_inner)
        
        # Compute the maximum truncation error
        te = h * f @ (np.array(c4) - np.array(c5))
        te_max = np.max(np.abs(te))
        
        # Compute the allowable truncation error
        ymax = np.max(np.abs(y))
        te_allowed = tolerance * max(ymax, 1.0)
        
        # Compute the fractional change in step size
        delta = (te_allowed / (te_max + np.finfo(float).eps)) ** (1/5)
        
        # If the truncation error is in bounds, then update the solution
        if te_max <= te_allowed:
            h = min(h, tf - t)
            t += h
            y = yi + h * f @ np.array(c5)
            tout.append(t)
            yout.append(y.copy())
        
        # Update the time step
        h = min(delta * h, 4 * h)
        if h < hmin:
            print(f'\n\n Warning: Step size fell below its minimum allowable value ({hmin}) at time {t}.\n\n')
            break

    return np.array(tout), np.array(yout)
