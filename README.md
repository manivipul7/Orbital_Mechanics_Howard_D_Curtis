# Orbital Mechanics – Study Scripts

This repository contains my personal study implementation and organization of scripts for classical orbital mechanics and spacecraft dynamics [Curtis, H.D., 2019. Orbital mechanics for engineering students. Butterworth-Heinemann.].  
The work is inspired by examples and exercises from “Orbital Mechanics for Engineering Students” by Howard D. Curtis; that book and its original materials remain the property of their respective copyright holder.

The code is a mix of MATLAB and Python and is organized roughly by chapter/topic.  
File names and example numbers (D.xx, Example_x_xx.m) follow the same labels used in the source material so it is easy to cross‑reference when studying.

> Note: This repository is intended for learning and experimentation only.  
> It is not an official redistribution of the textbook’s code or text.

---

## Chapter 1 – Dynamics of Point Masses

Time‑integration routines for point‑mass motion and basic ODE examples.

- D.1 – Introduction  
- D.2 – `rkf1_4.m`: Runge–Kutta time‑stepping (orders 1–4) for simple dynamics.  
  - `Example_1_18.m`: Sample problem using the RK routines.  
- D.3 – `heun.m`: Two‑stage predictor–corrector integrator (Heun method).  
  - `Example_1_19.m`: Example using Heun’s method on a test problem.  
- D.4 – `rk45.m`: Adaptive Runge–Kutta–Fehlberg 4(5) stepper for improved accuracy control.

---

## Chapter 2 – The Two‑Body Problem

Tools for propagating Keplerian orbits and exploring two‑body motion.

- `Example_1_20.m`: Additional worked example using the adaptive integrator.  
- D.5 – `twobody3d.m`: 3D two‑body propagation with sample initial conditions (Example 2.2 style).  
- D.6 – `orbit.m`: Relative two‑body motion between objects sharing a central body.  
- D.7 – `f_and_g_ta.m`: Lagrange f and g coefficients as a function of true‑anomaly change.  
  - `fDot_and_gDot_ta.m`: Time‑derivative forms of f and g for the same parameterization.  
- D.8 – `rv_from_r0v0_ta.m`: State‑vector update from an initial state and a change in true anomaly.  
  - `Example_2_13.m`: Worked case for the above routine.  
- D.9 – `bisect.m`: Generic bisection root‑finder used by several examples.  
  - `Example_2_16.m`: Demonstration of the bisection solver.  
- D.10 – `Example_2_18.m`: Translunar transfer modeled as a restricted three‑body trajectory.

---

## Chapter 3 – Orbital Position as a Function of Time

Implementations of Kepler equation solvers and universal‑variable propagation.

- D.11 – `kepler_E.m`: Newton–Raphson solver for Kepler’s equation in elliptic form.  
  - `Example_3_02.m`: Example using the elliptic Kepler solver.  
- D.12 – `kepler_H.m`: Hyperbolic Kepler equation solver using Newton’s method.  
  - `Example_3_05.m`: Hyperbolic‑orbit example.  
- D.13 – `stumpS.m`, `stumpC.m`: Stumpff functions used in the universal‑variable formulation.  
- D.14 – `kepler_U.m`: Universal‑variable Kepler solver covering multiple conic types.  
  - `Example_3_06.m`: Example of the universal‑variable approach.  
- D.15 – `f_and_g.m`: Lagrange f and g and their time derivatives using universal anomaly.  
- D.16 – `rv_from_r0v0.m`: State‑vector propagation from initial state and time of flight.  
  - `Example_3_07.m`: Sample propagation case.

---

## Chapter 4 – Orbits in Three Dimensions

Conversions between state vectors, orbital elements, and attitude‑related angles.

- D.17 – `ra_and_dec_from_r.m`: Right ascension and declination from a position vector.  
  - `Example_4_01.m`: Example of converting position to (RA, Dec).  
- D.18 – `coe_from_sv.m`: Classical orbital elements derived from a Cartesian state vector.  
  - `Example_4_03.m`: Example of state‑to‑elements conversion.  
- D.19 – `atan2d_0_360.m`: Angle utility returning results on the full [0°, 360°] range.  
- D.20 – `dcm_to_euler.m`: Euler‑angle sequence extracted from a direction cosine matrix (DCM).  
- D.21 – `dcm_to_ypr.m`: Yaw–pitch–roll angles from a DCM.  
- D.22 – `sv_from_coe.m`: State vector reconstruction from classical elements.  
  - `Example_4_07.m`: Example of elements‑to‑state conversion.  
- D.23 – `ground_track.m`: Basic ground‑track generation for a satellite given its elements.

---

## Chapter 5 – Preliminary Orbit Determination

Scripts for initial orbit estimation from position/angle data.

- D.24 – `gibbs.m`: Gibbs method implementation for three‑position orbit determination.  
  - `Example_5_01.m`: Sample orbit‑determination case using Gibbs.  
- D.25 – `lambert.m`: Lambert transfer solver for two‑point boundary value problems.  
  - `Example_5_02.m`: Example of Lambert transfer between two position vectors.  
- D.26 – `J0.m`: Julian day at 0 h UT, used as a time‑keeping utility.  
  - `Example_5_04.m`: Example of Julian day computation.  
- D.27 – `LST.m`: Local sidereal time evaluation for a given site and epoch.  
  - `Example_5_06.m`: Example of LST calculation.  
- D.28 – `rv_from_observe.m`: State estimation from range/angles and their time derivatives.  
  - `Example_5_10.m`: Sample tracking‑data orbit reconstruction.  
- D.29 – `gauss.m`: Gauss‑type preliminary orbit determination with refinement loop.  
  - `Example_5_11.m`: Example problem using the Gauss routine.

---

## Chapter 6 – Orbital Maneuvers

- D.30 – `integrate_thrust.m`: Finite‑burn propagation with constant thrust over a specified duration (Example 6.15‑type maneuver).

---

## Chapter 7 – Relative Motion and Rendezvous

- D.31 – `rva_relative.m`: Relative position, velocity, and acceleration in a comoving reference frame.  
  - `Example_7_01.m`: Example of relative motion between two spacecraft.  
- D.32 – `Example_7_02.m`: Plotting and visualizing relative trajectories.  
- D.33 – `Example_7_03.m`: Linearized relative‑motion solution with an elliptical reference orbit.

---

## Chapter 8 – Interplanetary Trajectories

- D.34 – `month_planet_names.m`: Utility to map numeric month/planet codes to their names.  
- D.35 – `planet_elements_and_sv.m`: Planet ephemeris approximation and heliocentric state vectors.  
  - `Example_8_07.m`: Example of computing a planet’s heliocentric state.  
- D.36 – `interplanetary.m`: Simple interplanetary transfer trajectory between two planets.  
  - `Example_8_08.m`: Example interplanetary mission profile.

---

## Chapter 9 – Lunar Trajectories

- D.37 – `simpsons_lunar_ephemeris.m`: Approximate lunar state as a function of time.  
- D.38 – `Example_9_03.m`: Numerical lunar‑trajectory example using the ephemeris.

---

## Chapter 10 – Introduction to Orbital Perturbations

- D.39 – `atmosphere.m`: Standard‑atmosphere model for basic drag and density estimates.  
- D.40 – `Example_10_01.m`: Orbit‑decay time under atmospheric drag (Cowell‑type approach).  
- D.41 – `Example_10_02.m`: J2‑perturbed orbit using an Encke‑style formulation.  
- D.42 – `Example_10_06.m`: Use of Gauss variational equations for J2‑driven element changes.  
- D.43 – `solar_position.m`: Geocentric solar position for a given epoch.  
- D.44 – `los.m`: Line‑of‑sight and eclipse check to determine if a satellite is in shadow.  
- D.45 – `Example_10_09.m`: Solar‑radiation‑pressure driven element evolution.  
- D.46 – `lunar_position.m`: Geocentric lunar position for a given epoch.  
  - `Example_10_10.m`: Example based on the lunar‑position routine.  
- D.47 – `Example_10_11.m`: Lunar‑gravity perturbation of orbital elements.  
- D.48 – `Example_10_12.m`: Solar‑gravity perturbation of orbital elements.

---

## Chapter 11 – Rigid Body Dynamics

- D.49 – `dcm_from_q.m`: Direction cosine matrix from a quaternion.  
- D.50 – `q_from_dcm.m`: Quaternion computed from a given DCM.  
- D.51 – `quat_rotate.m`: Quaternion‑based vector rotation utility.  
- D.52 – `Example_11_26.m`: Classic spinning‑top style rigid‑body example.

---

## Chapters 12–13 – Attitude and Rocket Dynamics

- **Chapter 12: Spacecraft Attitude Dynamics** – placeholder for additional attitude scripts.  
- **Chapter 13: Rocket Vehicle Dynamics**  
  - D.53 – `Example_13_03.m`: Gravity‑turn‑type ascent/trajectory example.

---

## Acknowledgements

- Curtis, H.D., 2019. Orbital mechanics for engineering students. Butterworth-Heinemann., which motivated the structure and many of the example problems in this study repository.
