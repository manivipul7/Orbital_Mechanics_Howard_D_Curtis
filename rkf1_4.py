import numpy as np

def rkf1_4(ode_function, tspan, y0, h, rk):
    # Determine which Runge-Kutta method is to be used:
    if rk == 1:
        n_stages = 1
        a = np.array([0])
        b = np.array([[0]])
        c = np.array([1])
    elif rk == 2:
        n_stages = 2
        a = np.array([0, 1])
        b = np.array([[0, 0], [1, 0]])
        c = np.array([1/2, 1/2])
    elif rk == 3:
        n_stages = 3
        a = np.array([0, 1/2, 1])
        b = np.array([[0, 0, 0], [1/2, 0, 0], [-1, 2, 0]])
        c = np.array([1/6, 2/3, 1/6])
    elif rk == 4:
        n_stages = 4
        a = np.array([0, 1/2, 1/2, 1])
        b = np.array([[0, 0, 0, 0], [1/2, 0, 0, 0], [0, 1/2, 0, 0], [0, 0, 1, 0]])
        c = np.array([1/6, 1/3, 1/3, 1/6])
    else:
        raise ValueError('The parameter rk must have the value 1, 2, 3 or 4.')

    t0, tf = tspan
    t = t0
    y = y0.astype(np.float64)  # Ensure y is float64
    tout = [t]
    yout = [y.copy()]
    
    f = np.zeros((len(y0), n_stages), dtype=np.float64)  # Ensure f is float64

    while t < tf:
        ti = t
        yi = y.copy()
        
        # Evaluate the time derivatives at the n_stages points within the current interval
        for i in range(n_stages):
            t_inner = ti + a[i] * h
            y_inner = yi.copy()
            for j in range(i):
                y_inner += h * b[i, j] * f[:, j]
            f[:, i] = ode_function(t_inner, y_inner)
        
        h = min(h, tf - t)
        t += h
        y = yi + h * np.dot(f, c)
        tout.append(t)
        yout.append(y.copy())

    return np.array(tout), np.array(yout)
