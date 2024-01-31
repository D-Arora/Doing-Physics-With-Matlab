# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 08:15:33 2024

@author: Owner
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  1 16:45:28 2022

@author: laurabarbara
"""

from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp as rk4

plt.rcParams["figure.figsize"] = (7,7)

#------------------------------------#
#                RK4                 #
#------------------------------------#

def fsp_derivatives(tau,y):
    ''' nondimentionalised eq. of motion written as coupled ODEs to be solved by RK4'''
    theta=y[0]
    omega=y[1]
    return [omega, -k*omega-np.sin(theta)+F*np.cos((1-eta)*tau)]

k=0.5
F=0.9
eta=1/3
tau=np.linspace(0,300,100000)

sol=rk4(fun=fsp_derivatives, t_span=[tau[0],tau[-1]], y0=[np.radians(0),0], t_eval=tau)

#------------------------------------#
#         wrapping the angle         #
#------------------------------------#

th_wrapped=np.arctan2(np.sin(sol.y[0]), np.cos(sol.y[0]))

#------------------------------------#
#             plotting               #
#------------------------------------#

plt.clf()

fig=plt.figure(1)
fig.suptitle('Non-linear periodic motion of FSP with ' + r'$F=$' + str(F) + ', ' + r'$k=$' + str(k) \
             + ', ' + r'$\eta=1/3$' + ' and ' + r'$~\theta(0)=0$', fontsize=16)

#------------------------------------#
#               fig 1                #
#------------------------------------#

f1=plt.subplot2grid((2,2), (0,0), colspan=2)
plt.scatter(sol.y[0][:], sol.y[1][:], c='b', s=0.5)
plt.grid()
f1.set_xlabel('Oscillation Amplitude, ' + r'$\theta$' + ' [rad]', fontsize=14)
f1.set_ylabel('Angular Velocity, ' + r'$\omega$' + ' [rad/s]', fontsize=14)
f1.set_xlim(-1*np.pi, np.pi)

#------------------------------------#
#               fig 2                #
#------------------------------------#

f2=plt.subplot2grid((2,2), (1,0), colspan=2)
plt.scatter(tau[:], th_wrapped[:], c='b', s=0.5)
plt.grid()
f2.set_ylabel('Oscillation Amplitude, ' + r'$\theta$' + ' [rad]', fontsize=14)
f2.set_xlabel(r'$\tau$', fontsize=14)
f2.set_xlim(tau[0], tau[-1])
f2.set_ylim(-1*np.pi, np.pi)

plt.tight_layout()
#plt.savefig('non linear.png', dpi=400)
plt.show()