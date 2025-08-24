# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 09:39:39 2025

@author: Owner
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# 1. Define the Hindmarsh-Rose model
def hindmarsh_rose(X, t, a, b, c, d, r, s, x0, I):
    x, y, z = X
    dxdt = y - a*x**3 + b*x**2 - z + I
    dydt = c - d*x**2 - y
    dzdt = r*(s*(x - x0) - z)
    return [dxdt, dydt, dzdt]

# Set fixed parameters
a, b, c, d, r, s, x0 = 1, 3, 1, 5, 0.006, 4, -1.6

# Define range for the control parameter (e.g., I)
I_values = np.linspace(1.0, 4.0, 100) # Example range

all_isis = []
for I in I_values:
    # Initial conditions
    X0 = [-1.5, -5.0, 0.0]
    t = np.linspace(0, 500, 10000) # Simulation time

    # Solve ODEs
    sol = odeint(hindmarsh_rose, X0, t, args=(a, b, c, d, r, s, x0, I))
    x_membrane = sol[:, 0]

    # 3. Spike detection and ISI calculation (simplified example)
    spike_times = []
    threshold = 0.0 # Example threshold for spike detection
    for i in range(1, len(x_membrane)):
        if x_membrane[i] > threshold and x_membrane[i-1] <= threshold:
            spike_times.append(t[i])

    isis = [spike_times[i] - spike_times[i-1] for i in range(1, len(spike_times))]
    
    # Store ISIs with corresponding I value
    for isi in isis:
        all_isis.append((I, isi))

# 5. Plotting the bifurcation diagram
I_plot = [item[0] for item in all_isis]
isi_plot = [item[1] for item in all_isis]

plt.figure(figsize=(10, 6))
plt.scatter(I_plot, isi_plot, s=1, alpha=0.5)
plt.xlabel("External Current (I)")
plt.ylabel("Interspike Interval (ISI)")
plt.title("ISI Bifurcation Diagram of Hindmarsh-Rose Model")
plt.grid(True)
plt.show()