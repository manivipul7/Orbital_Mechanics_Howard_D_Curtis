# Orbital Mechanics

Python code based on the programs provided in "Orbital Mechanics for Engineering Students" by Howard D. Curtis.

| Appendix | File Name                   | Description |
|----------|-----------------------------|-------------|
| **Chapter 1: Dynamics of Point Masses**| |
| D.1      | -                           | Introduction |
| D.2      | rkf1_4.m                    | Algorithm 1.1: Numerical integration using Runge-Kutta methods RK1, RK2, RK3, or RK4. |
|          | Example_1_18.m              | Example of Algorithm 1.1. |
| D.3      | heun.m                      | Algorithm 1.2: Numerical integration using Heun’s method. |
|          | Example_1_19.m              | Example of Algorithm 1.2. |
| D.4      | rk45.m                      | Algorithm 1.3: Numerical integration using the Runge-Kutta-Fehlberg 4(5) method with adaptive step size control. |
| **Chapter 2: The Two-body Problem** ||
|          | Example_1_20.m              | Example of Algorithm 1.3. |
| D.5      | twobody3d.m                 | Algorithm 2.1: Numerical solution for two-body motion in 3D. Includes data for Example 2.2. |
| D.6      | orbit.m                     | Algorithm 2.2: Numerical solution for two-body relative motion. Includes data for Example 2.3. |
| D.7      | f_and_g_ta.m                | Calculates Lagrange coefficients f and g in terms of change in true anomaly. |
|          | fDot_and_gDot_ta.m          | Calculates derivatives of Lagrange coefficients f and g in terms of change in true anomaly. |
| D.8      | rv_from_r0v0_ta.m           | Algorithm 2.3: Calculate state vector given initial state vector and change in true anomaly. |
|          | Example_2_13.m              | Example of Algorithm 2.3. |
| D.9      | bisect.m                    | Algorithm 2.4: Find root of a function using bisection method. |
|          | Example_2_16.m              | Example of Algorithm 2.4. |
| D.10     | Example_2_18.m              | Translunar trajectory as a circular restricted three-body problem. |
| **Chapter 3: Orbital Position as a Function of Time** ||
| D.11     | kepler_E.m                  | Algorithm 3.1: Solution of Kepler’s equation by Newton’s method. |
|          | Example_3_02.m              | Example of Algorithm 3.1. |
| D.12     | kepler_H.m                  | Algorithm 3.2: Solution of Kepler’s equation for hyperbola using Newton’s method. |
|          | Example_3_05.m              | Example of Algorithm 3.2. |
| D.13     | stumpS.m                    | Calculation of the Stumpff function S(z) and C(z). |
|          | stumpC.m                    | Calculation of the Stumpff function C(z). |
| D.14     | kepler_U.m                  | Algorithm 3.3: Solution of the universal Kepler’s equation using Newton’s method. |
|          | Example_3_06.m              | Example of Algorithm 3.3. |
| D.15     | f_and_g.m                   | Calculation of Lagrange coefficients f and g and their time derivatives in terms of change in universal anomaly. |
| D.16     | rv_from_r0v0.m              | Algorithm 3.4: Calculation of the state vector given the initial state vector and the time lapse Δt. |
|          | Example_3_07.m              | Example of Algorithm 3.4. |
| **Chapter 4: Orbits in Three Dimensions** ||
| D.17     | ra_and_dec_from_r.m         | Algorithm 4.1: Obtain right ascension and declination from the position vector. |
|          | Example_4_01.m              | Example of Algorithm 4.1. |
| D.18     | coe_from_sv.m               | Algorithm 4.2: Calculation of the orbital elements from the state vector. |
|          | Example_4_03.m              | Example of Algorithm 4.2. |
| D.19     | atan2d_0_360.m              | Calculation of atan2(y/x) in the range 0° to 360°. (MATLAB’s atan2d result lies in the range 0° to 180°.) |
| D.20     | dcm_to_euler.m              | Algorithm 4.3: Obtain the classical Euler angle sequence from a Direction Cosine Matrix (DCM). |
| D.21     | dcm_to_ypr.m                | Algorithm 4.4: Obtain yaw, pitch, and roll angles from a DCM. |
| D.22     | sv_from_coe.m               | Algorithm 4.5: Calculation of the state vector from the orbital elements. |
|          | Example_4_07.m              | Example of Algorithm 4.5. |
| D.23     | ground_track.m              | Algorithm 4.6: Calculate the ground track of a satellite from its orbital elements. Contains data for Example 4.12. |
| **Chapter 5: Preliminary Orbit Determination** ||
| D.24     | gibbs.m                     | Algorithm 5.1: Gibbs’ method of preliminary orbit determination. |
|          | Example_5_01.m              | Example of Algorithm 5.1. |
| D.25     | lambert.m                   | Algorithm 5.2: Solution of Lambert’s problem. |
|          | Example_5_02.m              | Example of Algorithm 5.2. |
| D.26     | J0.m                        | Calculation of Julian day number at 0 hr UT. |
|          | Example_5_04.m              | Example of Julian day calculation. |
| D.27     | LST.m                       | Algorithm 5.3: Calculation of Local Sidereal Time (LST). |
|          | Example_5_06.m              | Example of Algorithm 5.3. |
| D.28     | rv_from_observe.m           | Algorithm 5.4: Calculation of the state vector from measurements of range, angular position, and rates. |
|          | Example_5_10.m              | Example of Algorithm 5.4. |
| D.29     | gauss.m                     | Algorithms 5.5 and 5.6: Gauss’ method of preliminary orbit determination with iterative improvement. |
|          | Example_5_11.m              | Example of Algorithms 5.5 and 5.6. |
| **Chapter 6: Orbital Maneuvers** ||
| D.30     | integrate_thrust.m          | Calculate the state vector at the end of a finite time, constant thrust delta-v maneuver. Contains data for Example 6.15. |
| **Chapter 7: Relative Motion and Rendezvous** ||
| D.31     | rva_relative.m              | Algorithm 7.1: Find the position, velocity, and acceleration of B relative to A’s comoving frame. |
|          | Example_7_01.m              | Example of Algorithm 7.1. |
| D.32     | Example_7_02.m              | Plot the position of one spacecraft relative to another. |
| D.33     | Example_7_03.m              | Solve the linearized equations of relative motion with an elliptical reference orbit. |
| **Chapter 8: Interplanetary Trajectories** ||
| D.34     | month_planet_names.m        | Convert the numerical designation of a month or a planet into its name. |
| D.35     | planet_elements_and_sv.m    | Algorithm 8.1: Calculation of the heliocentric state vector of a planet at a given epoch. |
|          | Example_8_07.m              | Example of Algorithm 8.1. |
| D.36     | interplanetary.m            | Algorithm 8.2: Calculate the spacecraft trajectory from planet 1 to planet 2. |
|          | Example_8_08.m              | Example of Algorithm 8.2. |
| **Chapter 9: Lunar Trajectories** ||
| D.37     | simpsons_lunar_ephemeris.m  | Calculation of lunar state vector vs. time. |
| D.38     | Example_9_03.m              | Numerical calculation of lunar trajectory. |
| **Chapter 10: Introduction to Orbital Perturbations** ||
| D.39     | atmosphere.m                | US Standard Atmosphere 1976. |
| D.40     | Example_10_01.m             | Time for orbit decay using Cowell’s method. |
| D.41     | Example_10_02.m             | J2 perturbation of an orbit using Encke’s method. |
| D.42     | Example_10_06.m             | Using Gauss’ variational equations to assess the J2 effect on orbital elements. |
| D.43     | solar_position.m            | Algorithm 10.2: Calculate the geocentric position of the sun at a given epoch. |
| D.44     | los.m                       | Algorithm 10.3: Determine whether or not a satellite is in Earth’s shadow. |
| D.45     | Example_10_09.m             | Use Gauss’ variational equations to determine the effect of solar radiation pressure on Earth satellite’s orbital parameters. |
| D.46     | lunar_position.m            | Algorithm 10.4: Calculate the geocentric position of the moon at a given epoch. |
|          | Example_10_10.m             | Example of Algorithm 10.4. |
| D.47     | Example_10_11.m             | Use Gauss’ variational equations to determine the effect of lunar gravity on Earth satellite’s orbital parameters. |
| D.48     |  Example_10_12.m            |  Use the Gauss’ variational equations to determine the effect of solar gravity on an earth satellite’s orbital parameters. |
| **Chapter 11: Rigid Body Dynamics** ||
| D.49     |  dcm_from_q.m               | Algorithm 11.1: Calculate the direction cosine matrix from the quaternion. |
| D.50     |  q_from_dcm.m               | Algorithm 11.2: Calculate the quaternion from the direction cosine matrix. |
| D.51     |   quat_rotate.m             | Quaternion vector rotation operation (Eq. 11.160). |
| D.52     |   Example_11_26.m           | Solution of the spinning top problem. |
| **Chapter 12: Spacecraft Attitude Dynamics** ||
| **Chapter 13: Rocket Vehicle Dynamics** ||
| D.53     |   Example_13_03.m           | Example 11.3: Calculation of a gravity turn trajectory. |







