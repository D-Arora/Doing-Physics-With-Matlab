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
#from matplotlib.animation import FuncAnimation, PillowWriter 
from scipy.integrate import solve_ivp
import time

from matplotlib.animation import FuncAnimation

from IPython.display import HTML



tStart = time.time()
plt.close('all')

def lorenz(t, state, sigma=10, beta=8/3, rho=28):
    x, y, z = state
    dxdt = sigma * (y - x)
    dydt = x * (rho - z) - y
    dzdt = x * y - beta * z
    return [dxdt, dydt, dzdt]


# Multiple close initial conditions
initial_conditions = [
    [5., 5., 5.],
    [5.1, 5.1, 5.2],
    [4.8, 5., 4.8]
]

# Longer simulation
t_span = (0, 60)
t_eval = np.linspace(*t_span, 3000)

def integrate_lorenz(initial_state, t_span, t_eval):
    sol = solve_ivp(lorenz, t_span, initial_state, t_eval=t_eval, method='RK45')
    return sol.y

trajectories = [integrate_lorenz(init, t_span, t_eval) for init in initial_conditions]



# Create 3D plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(projection='3d')
ax.set_xlim([-25, 25])
ax.set_ylim([-35, 35])
ax.set_zlim([0, 50])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Lorenz Attractor: Animation of 3 Initial Conditions')

colors = ['blue', 'red', 'green']
lines = [ax.plot([], [], [], lw=1, color=colors[i], markersize=0.5)[0] for i in range(3)]

# Animation update function
def update(num):
    for i in range(3):
        lines[i].set_data(trajectories[i][0][:num], trajectories[i][1][:num])
        lines[i].set_3d_properties(trajectories[i][2][:num])
    return lines

# Run animation
anim = FuncAnimation(fig, update, frames=1000, interval=30, blit=True)

# Display
HTML(anim.to_jshtml())
#anim.save('agLE.gif', fps = 5)

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


