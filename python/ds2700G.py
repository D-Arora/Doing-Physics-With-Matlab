# -*- coding: utf-8 -*-
'''
ds2700A.py      Nov 2025
DYNAMICAL SYSTEMS
ANIMATION: LORENZ EQUATIONS

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDS.ds2700A.pdf
'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, ones, sqrt, real, imag 
from scipy.integrate import odeint
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

tStart = time.time()
plt.close('all')

#%%
def lorenz_deriv(state, t):
    x, y, z = state
    dx_dt = s * (y - x)
    dy_dt = x * (r - z) - y
    dz_dt = x * y - b * z
    return [dx_dt, dy_dt, dz_dt]
    
# Model parameters, Initial conditions
r = 28
s = 10
b = 8/3
f = 3
initial_state = [2, 2, 0]

# Time points
N = 1988
t = np.linspace(0, 40, N) # Adjust time range and number of points as needed

# Solve the system
solution = odeint(lorenz_deriv, initial_state, t)
x, y, z = solution[:, 0], solution[:, 1], solution[:, 2] 

initial_state = [-2, -2, 0] 

solution = odeint(lorenz_deriv, initial_state, t)
X, Y, Z = solution[:, 0], solution[:, 1], solution[:, 2]   

# Fixed points other than the Origin
xEp,xEm,yEp,yEm,zE = 0,0,0,0,0
if r > 1:
   ev1 = zeros(N); ev2 = ev1; ev3 = ev1 
   xEp = sqrt(b*(r-1)); xEm = - xEp
   yEp = xEp; yEm = xEm; zE = r-1


#%% Animated plot
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,5)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x'); ax.set_ylabel('y');ax.set_zlabel('z')
ax.set_title('r = %0.3f' %r,fontsize = 12)
line, =  ax.plot([], [], [], 'b', lw=0.5) # 'b' for blue, lw for linewidth
line1, = ax.plot([], [], [], 'r', lw=0.5)
line2, = ax.plot([], [], [], 'b', lw=4)
line3, = ax.plot([], [], [], 'r', lw=4)
ax.plot3D(xEp,yEp,zE,'ko',ms = 7)
ax.plot3D(xEm,yEm,zE,'ko',ms = 7)
ax.plot3D(0,0,0,'ko',ms = 7)

# Set initial plot limits (optional, can be adjusted for better view)
ax.set_xlim([-25, 25])
ax.set_ylim([-35, 35])
ax.set_zlim([0, 50])
ax.view_init(elev=23, azim=-50, roll = -7)


def init():
    line.set_data([], [])
    line.set_3d_properties([])
    line1.set_data([], [])
    line1.set_3d_properties([])
    line2.set_data([], [])
    line2.set_3d_properties([])
    line3.set_data([], [])
    line3.set_3d_properties([])
    ax.plot3D(0,0,0,'ko',ms = 7)
    return line, line1, line2, line3,

def update(frame):
# Update the data for the line plot for each frame
    line.set_data(x[:frame], y[:frame])
    line.set_3d_properties(z[:frame])
    line1.set_data(X[:frame], Y[:frame])
    line1.set_3d_properties(Z[:frame])
    line2.set_data(x[frame-f:frame], y[frame-f:frame])
    line2.set_3d_properties(z[frame-f:frame]) 
    line3.set_data(X[frame-f:frame], Y[frame-f:frame])
    line3.set_3d_properties(Z[frame-f:frame]) 
    ax.plot3D(0,0,0,'ko',ms = 7)
    time.sleep(0.0001)
    return line, line1, line2, line3,

# Create the animation

anim = FuncAnimation(fig, update, frames=len(t), init_func=init,
                     interval=30, blit=True, repeat = False)

anim.save('agLE.gif', writer = 'pillow', fps = 120)



#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


