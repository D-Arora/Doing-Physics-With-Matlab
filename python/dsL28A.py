# -*- coding: utf-8 -*-
"""
dsL28A.py
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import time

#%%
tStart = time.time()
plt.close('all')




# Parameters
g = 9.81
L = 1.0
gamma = 0.2          # damping
A = 1.2              # drive amplitude
omega = 2.0/3.0      # drive frequency

def pendulum(t, y):
    theta, omega_theta = y
    dtheta = omega_theta
    domega = -(g/L)*np.sin(theta) - gamma*omega_theta + A*np.cos(omega*t)
    return [dtheta, domega]

# Integration settings
T_drive = 2*np.pi/omega          # drive period
n_periods = 2000                 # simulate many periods
t_span = (0, n_periods*T_drive)
t_eval = np.linspace(*t_span, n_periods*200)  # fine sampling

# Integrate
y0 = [0.2, 0.0]   # initial angle and angular velocity
sol = solve_ivp(pendulum, t_span, y0, t_eval=t_eval)
t = sol.t
theta = sol.y[0]
omega_theta = sol.y[1]

# Sample once per period (stroboscopic map)
k = np.arange(200, n_periods)      # skip initial transients
t_sample = k * T_drive

# Interpolate solution at sampling times
theta_sample = np.interp(t_sample, t, theta)
omega_sample = np.interp(t_sample, t, omega_theta)

# Wrap angle into [-pi, pi] for nicer plot
theta_sample = (theta_sample + np.pi) % (2*np.pi) - np.pi

plt.figure(figsize=(5, 5))
plt.scatter(theta_sample, omega_sample, s=1, color='k')
plt.xlabel(r'$\theta$')
plt.ylabel(r'$\dot{\theta}$')
plt.title('Poincaré section of driven pendulum')
plt.grid(True)
plt.tight_layout()
plt.show()

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(np.round(tExe,2))