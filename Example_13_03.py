import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def example_13_03():
    deg = np.pi / 180  # Convert degrees to radians
    g0 = 9.81  # Sea-level acceleration of gravity (m/s^2)
    Re = 6378e3  # Radius of the earth (m)
    hscale = 7.5e3  # Density scale height (m)
    rho0 = 1.225  # Sea level density of atmosphere (kg/m^3)
    
    diam = 196.85 / 12 * 0.3048  # Vehicle diameter (m)
    A = np.pi / 4 * (diam) ** 2  # Frontal area (m^2)
    CD = 0.5  # Drag coefficient (assumed constant)
    m0 = 149912 * 0.4536  # Lift-off mass (kg)
    n = 7  # Mass ratio
    T2W = 1.4  # Thrust to weight ratio
    Isp = 390  # Specific impulse (s)
    
    mfinal = m0 / n  # Burnout mass (kg)
    Thrust = T2W * m0 * g0  # Rocket thrust (N)
    m_dot = Thrust / Isp / g0  # Propellant mass flow rate (kg/s)
    mprop = m0 - mfinal  # Propellant mass (kg)
    tburn = mprop / m_dot  # Burn time (s)
    hturn = 130  # Height at which pitchover begins (m)
    
    t0 = 0  # Initial time for the numerical integration
    tf = tburn  # Final time for the numerical integration
    tspan = (t0, tf)  # Range of integration
    
    # Initial conditions:
    v0 = 0  # Initial velocity (m/s)
    gamma0 = 89.85 * deg  # Initial flight path angle (rad)
    x0 = 0  # Initial downrange distance (km)
    h0 = 0  # Initial altitude (km)
    vD0 = 0  # Initial value of velocity loss due to drag (m/s)
    vG0 = 0  # Initial value of velocity loss due to gravity (m/s)
    
    # Initial conditions vector:
    f0 = [v0, gamma0, x0, h0, vD0, vG0]
    
    # Solve the system of equations df/dt = f(t) using Runge-Kutta method
    sol = solve_ivp(rates, tspan, f0, args=(m0, m_dot, tburn, Thrust, g0, Re, hscale, rho0, A, CD, hturn), 
                method='RK45', t_eval=np.linspace(t0, tf, 1000), rtol=1e-10, atol=1e-14)
    
    t = sol.t
    f = sol.y.T
    
    # Extract the solution components
    v = f[:, 0] * 1.e-3  # Velocity (km/s)
    gamma = f[:, 1] / deg  # Flight path angle (degrees)
    x = f[:, 2] * 1.e-3  # Downrange distance (km)
    h = f[:, 3] * 1.e-3  # Altitude (km)
    vD = -f[:, 4] * 1.e-3  # Velocity loss due to drag (km/s)
    vG = -f[:, 5] * 1.e-3  # Velocity loss due to gravity (km/s)
    
    # Dynamic pressure vs time:
    q = []
    M = []
    for i in range(len(t)):
        Rho = rho0 * np.exp(-h[i] * 1000 / hscale)  # Air density (kg/m^3)
        q.append(0.5 * Rho * (v[i] * 1.e3) ** 2)  # Dynamic pressure (Pa)
        a = 340.3  # Approximation of speed of sound (m/s) at sea level, should ideally use atmosphere model
        M.append(1000 * v[i] / a)  # Mach number
    
    q = np.array(q)
    M = np.array(M)
    
    # Maximum dynamic pressure and corresponding time, speed, altitude, and Mach number:
    maxQ = np.max(q)
    imax = np.argmax(q)
    tQ = t[imax]
    vQ = v[imax]
    hQ = h[imax]
    aQ = 340.3  # Approximation, replace with atmospheric model if available
    MQ = 1000 * vQ / aQ
    
    output(gamma0, deg, hturn, tburn, maxQ, tQ, vQ, hQ, MQ, v, gamma, h, x, vD, vG)
    
    # Plotting results
    plt.figure('Trajectory and Dynamic Pressure')
    plt.subplot(2, 1, 1)
    plt.plot(x, h)
    plt.title('(a) Altitude vs Downrange Distance')
    plt.axis('equal')
    plt.xlabel('Downrange Distance (km)')
    plt.ylabel('Altitude (km)')
    plt.grid()
    
    plt.subplot(2, 1, 2)
    plt.plot(h, q * 9.869e-6)
    plt.title('(b) Dynamic Pressure vs Altitude')
    plt.xlabel('Altitude (km)')
    plt.ylabel('Dynamic pressure (atm)')
    plt.grid()
    
    plt.show()

def rates(t, y, m0, m_dot, tburn, Thrust, g0, Re, hscale, rho0, A, CD, hturn):
    # Extract variables
    v, gamma, x, h, vD, vG = y
    
    # When time t exceeds the burn time, set the thrust and the mass flow rate equal to zero
    if t < tburn:
        m = m0 - m_dot * t  # Current vehicle mass
        T = Thrust  # Current thrust
    else:
        m = m0 - m_dot * tburn  # Current vehicle mass
        T = 0  # Current thrust
    
    g = g0 / (1 + h / Re) ** 2  # Gravitational variation with altitude h
    rho = rho0 * np.exp(-h / hscale)  # Exponential density variation with altitude
    D = 0.5 * rho * v ** 2 * A * CD  # Drag
    
    # Define the first derivatives of v, gamma, x, h, vD, and vG
    if h <= hturn:
        gamma_dot = 0
        v_dot = T / m - D / m - g
        x_dot = 0
        h_dot = v
        vG_dot = -g
    else:
        v_dot = T / m - D / m - g * np.sin(gamma)
        gamma_dot = -1 / v * (g - v ** 2 / (Re + h)) * np.cos(gamma)
        x_dot = Re / (Re + h) * v * np.cos(gamma)
        h_dot = v * np.sin(gamma)
        vG_dot = -g * np.sin(gamma)
    
    vD_dot = -D / m
    
    # Return derivatives as a list
    return [v_dot, gamma_dot, x_dot, h_dot, vD_dot, vG_dot]

def output(gamma0, deg, hturn, tburn, maxQ, tQ, vQ, hQ, MQ, v, gamma, h, x, vD, vG):
    print('\n\n -----------------------------------\n')
    print(f' Initial flight path angle = {gamma0 / deg:10.3f} deg')
    print(f' Pitchover altitude        = {hturn:10.3f} m')
    print(f' Burn time                 = {tburn:10.3f} s')
    print(f' Maximum dynamic pressure  = {maxQ * 9.869e-6:10.3f} atm')
    print(f'    Time                   = {tQ / 60:10.3f} min')
    print(f'    Speed                  = {vQ:10.3f} km/s')
    print(f'    Altitude               = {hQ:10.3f} km')
    print(f'    Mach Number            = {MQ:10.3f}')
    print(' At burnout:')
    print(f'    Speed                  = {v[-1]:10.3f} km/s')
    print(f'    Flight path angle      = {gamma[-1]:10.3f} deg')
    print(f'    Altitude               = {h[-1]:10.3f} km')
    print(f'    Downrange distance     = {x[-1]:10.3f} km')
    print(f'    Drag loss              = {vD[-1]:10.3f} km/s')
    print(f'    Gravity loss           = {vG[-1]:10.3f} km/s')
    print('\n -----------------------------------\n')

if __name__ == '__main__':
    example_13_03()
