# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 21:33:51 2025

@author: Owner
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


plt.close('all')
# 1. Define the SEIR model differential equations
def seir_model(y, t, N, beta, alpha, gamma):
    S, E, I, R = y
    dSdt = -beta * S * I / N
    dEdt = beta * S * I / N - alpha * E
    dIdt = alpha * E - gamma * I
    dRdt = gamma * I
    return dSdt, dEdt, dIdt, dRdt

#plt.close('all')

# 2. Set initial parameters and conditions
N = 1000          # Total population
I0 = 1            # Initial number of infected individuals
E0 = 0            # Initial number of exposed individuals
R0 = 0            # Initial number of recovered individuals
S0 = N - I0 - E0 - R0 # Initial number of susceptible individuals

# Model parameters (example for a specific disease like COVID-19 with social distancing u=0.1)
t_incubation = 5.1 # days
t_infective = 3.3  # days
alpha = 1.0 / t_incubation # Rate from exposed to infectious
gamma = 1.0 / t_infective  # Recovery rate
R_naught = 2.4     # Basic reproduction number (R0, R naught)
beta = R_naught * gamma # Transmission rate

# Time vector (days)
t = np.linspace(0, 100, 100)

# Initial conditions vector
y0 = S0, E0, I0, R0

# 3. Integrate the SEIR equations over the time grid, t
result = odeint(seir_model, y0, t, args=(N, beta, alpha, gamma))
S, E, I, R = result.T

N = S+E+I+R
# 4. Plot the results
plt.figure(figsize=(10, 6))
plt.plot(t, S, 'b', label='Susceptible')
plt.plot(t, E, 'y', label='Exposed')
plt.plot(t, I, 'r', label='Infected')
plt.plot(t, R, 'g', label='Recovered')
plt.plot(t,N, 'k')
plt.xlabel('Time (days)')
plt.ylabel('Number of individuals')
plt.title('SEIR Model Simulation')
plt.legend()
plt.grid(True)
plt.show()