# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 16:47:27 2025


"""

import numpy as np

def calculate_mle(system_function, initial_state, d0, num_steps, dt, rescale_interval):
        # system_function: function defining the dynamics (e.g., Lorenz system)
        # initial_state: initial condition for the reference trajectory
        # d0: initial separation distance
        # num_steps: total number of integration steps
        # dt: time step
        # rescale_interval: steps between rescaling

        u1 = np.array(initial_state)
        # Initialize a perturbed state slightly offset from u1
        u2 = u1 + d0 * np.random.rand(len(u1)) # Example: random perturbation

        lyapunov_sum = 0
        for i in range(num_steps):
            # Evolve both trajectories
            u1_new = system_function(u1, dt)
            u2_new = system_function(u2, dt)

            if (i + 1) % rescale_interval == 0:
                # Calculate new distance
                d_new = np.linalg.norm(u2_new - u1_new)
                # Rescale u2_new to maintain d0 distance from u1_new
                u2_new = u1_new + d0 * (u2_new - u1_new) / d_new
                lyapunov_sum += np.log(d_new / d0)

            u1 = u1_new
            u2 = u2_new

        mle = lyapunov_sum / (num_steps * dt)
        return mle

    # Example: A placeholder for a system_function (e.g., for Lorenz system)
def lorenz_system(state, dt, sigma=10, rho=28, beta=8/3):
        x, y, z = state
        dx = sigma * (y - x)
        dy = x * (rho - z) - y
        dz = x * y - beta * z
        return np.array([x + dx * dt, y + dy * dt, z + dz * dt])

    # Example usage:
initial_state = [0.1, 0.0, 0.0]
mle_estimate = calculate_mle(lorenz_system, initial_state, d0=1e-9, num_steps=100000, dt=0.01, rescale_interval=100)

print(f"Estimated Maximum Lyapunov Exponent: {mle_estimate}")