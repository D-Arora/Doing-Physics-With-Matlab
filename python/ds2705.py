# -*- coding: utf-8 -*-
'''
ds2501.py      Oct 2025
DYNAMICAL SYSTEMS
PITCHFORK BIFURCATIONS IN PLANR SYSTEMS

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDS.ds2501.pdf
'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, ones, sqrt, real, imag 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import find_peaks
from numpy.linalg import eig
import time

tStart = time.time()
plt.close('all')

def lorenz(state, t, sigma, rho, beta):
    x, y, z = state
    dxdt = sigma * (y - x)
    dydt = x * (rho - z) - y
    dzdt = x * y - beta * z
    return [dxdt, dydt, dzdt]

def lorenz_jacobian(state, sigma, rho, beta):
    x, y, z = state
    J = np.array([[-sigma, sigma, 0], [rho - z, -1, -x], [y, x, -beta]])
    return J 


def calculate_lyapunov_exponent(sigma, rho, beta, initial_state, dt, num_steps, perturbation_magnitude=1e-9):
        n = len(initial_state)
        # Initial state and a set of orthonormal perturbation vectors
        state = np.array(initial_state)
        perturbations = np.identity(n) * perturbation_magnitude

        lyapunov_exponents = np.zeros(n)
        time = 0.0

        for _ in range(num_steps):
            # Integrate the main trajectory
            sol_main = odeint(lorenz, state, [0, dt], args=(sigma, rho, beta))
            state = sol_main[1]

            # Integrate the perturbations using the Jacobian
            J = lorenz_jacobian(state, sigma, rho, beta)
            for i in range(n):
                # Approximate perturbation evolution: d(delta_x)/dt = J * delta_x
                perturbations[:, i] = perturbations[:, i] + dt * J @ perturbations[:, i]

            # Gram-Schmidt orthonormalization
            for i in range(n):
                # Subtract projections onto previous vectors
                for j in range(i):
                    perturbations[:, i] -= np.dot(perturbations[:, i], perturbations[:, j]) * perturbations[:, j]
                # Normalize and record growth
                norm = np.linalg.norm(perturbations[:, i])
                lyapunov_exponents[i] += np.log(norm / perturbation_magnitude) # Adjust for initial magnitude
                perturbations[:, i] = perturbations[:, i] / norm * perturbation_magnitude # Renormalize

            time += dt

        return lyapunov_exponents / time

sigma = 10.0
rho = 28.0
beta = 8/3.0
initial_state = [0.1, 0.0, 0.0]
dt = 0.01
num_steps = 100000

lyapunov_spectrum = calculate_lyapunov_exponent(sigma, rho, beta, initial_state, dt, num_steps)
print("Lyapunov Exponents:", lyapunov_spectrum)


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


