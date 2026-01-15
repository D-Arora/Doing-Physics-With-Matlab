# -*- coding: utf-8 -*-
"""
dsL28B.py
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

import time

#%%
tStart = time.time()
plt.close('all')


# 1. Define the dynamical system equations
def system_ode(t, X, *args):
    # X = [x, x_dot, ...]
    # Return derivatives of X
    pass

# 2. Define a function or method to detect the Poincaré section crossing
def poincare_section_event(t, X, *args):
    # Return 0 when the condition is met (e.g., x = 0 or a specific phase)
    # The solver uses this to find the exact time of intersection
    pass
poincare_section_event.terminal = False # Don't stop integration
poincare_section_event.direction = 1  # Detect crossing in specific direction (e.g., upward)

# 3. Integrate and collect intersection points
initial_conditions = [1]
T_max = 10
t_span = (0, T_max)
# Use solve_ivp with the events parameter to get precise intersection times
solution = solve_ivp(system_ode, t_span, initial_conditions, events=poincare_section_event, dense_output=True)

# 4. Extract and plot the points
intersection_points = solution.y_events[0] # Array of state variables at each event
plt.scatter(intersection_points[:, 0], intersection_points[:, 1], s=5)
plt.xlabel('X at section')
plt.ylabel('Y at section')
plt.show()

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(np.round(tExe,2))
