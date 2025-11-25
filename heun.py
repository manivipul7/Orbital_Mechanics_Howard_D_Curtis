import numpy as np

def heun(ode_function, tspan, y0, h):
    tol = 1.e-6
    itermax = 100
    t0, tf = tspan
    t = t0
    y = y0.astype(np.float64)  # Ensure y is float64
    tout = [t]
    yout = [y.copy()]

    while t < tf:
        h = min(h, tf - t)
        t1 = t
        y1 = y.copy()
        f1 = ode_function(t1, y1)
        y2 = y1 + f1 * h
        t2 = t1 + h
        err = tol + 1
        iter = 0

        while err > tol and iter <= itermax:
            y2p = y2.copy()
            f2 = ode_function(t2, y2p)
            favg = (f1 + f2) / 2
            y2 = y1 + favg * h
            err = np.max(np.abs((y2 - y2p) / (y2 + np.finfo(float).eps)))
            iter += 1

        if iter > itermax:
            print(f'\n Maximum no. of iterations ({itermax}) exceeded at time = {t}')
            print(' in function \'heun.\'\n')
            break

        t += h
        y = y2
        tout.append(t)
        yout.append(y.copy())

    return np.array(tout), np.array(yout)
