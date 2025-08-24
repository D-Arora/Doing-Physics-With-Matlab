# -*- coding: utf-8 -*-
'''
bp002.py      July 2025

COMPLEX SYSTEMS
   Neuron Models with Adapting Feedback Synapses


Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
  https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/bs001.pdf

'''

#%% Libraries
import numpy as np
from numpy import pi, sqrt, linspace, zeros, array, real, imag, exp, sin 
from numpy.linalg import eig
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import random

tStart = time.time()

#plt.close('all')

#%% FUNCTIONS  Solve ODE for x,y    
def lorenz(t, state): 
    x, y = state
    # u = a*x; s = b*y
    # phiU = 3*u*exp(-u**2/2)
    # phiS = 3*s*exp(-s**2/2)
    dx = y - a*x**3 + b*x**2 + I0
    dy = c - d*x**2 - y
    return [dx, dy]  

# def fn_phi(u):
#     fn = 3*u*exp(-u**2/2)
#     return fn 
            
#%% INPUTS
a = 1
b = 3
c = 1
d = 5
I0 = 0
w = 2*pi

u0 = [0,1]
N = 9999
xS = zeros(N)
yS = zeros(N)
t = linspace(0,150,N)
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1]

plt.plot(xS[-520:-1],yS[-520:-1],'bo')
plt.plot(xS,yS)
plt.plot(t,xS,'k')
plt.grid()


#%% Console output
# print(' ')
# print('Critical point 0(xC, yC)')
# print('   (%0.0f  ' % xC[0] + ', %0.3f)' %yC[0] ) 
# print('Jacobian matrix J0')
# print(J0)
# print('Eigenvalues J0' )
# print('   (%0.3f  ' % V0[0] + ', %0.3f)' %V0[1] ) 
# print('Eigenfunctions J0')
# print(F0)

# print('Critical point 1(xC, yC)')
# print('   (%0.0f  ' % xC[1] + ', %0.3f)' %yC[1] ) 
# print('Jacobian matrix J1')
# print(J1)
# print('Eigenvalues J1' )
# print('   (%0.3f  ' % V1[0] + ', %0.3f)' %V1[1] ) 
# print('Eigenfunctions J1')
# print(F1)



#%% GRAPHICS
# FIGURE 1: x vs xDot   y vs yDot
# plt.rcParams['font.size'] = 12
# plt.rcParams["figure.figsize"] = (7,3)
# fig1, ax = plt.subplots(nrows=1, ncols=3)

# a = array([0.5, 1, 2])
# xP = u
# for C in range(3):
#     yP = fn_phi(a[C]*u)  
#     ax[C].set_xlabel('u',color= 'black',fontsize = 12)
#     ax[C].set_ylabel('$\phi(u)$',color = 'black',fontsize = 12)
#     ax[C].set_title('a = %0.2f' % a[C],color = 'black',fontsize = 12)
#     ax[C].set_xlim([-10,10])
#     ax[C].plot(xP,yP,'b', lw = 2)
#     ax[C].grid()
# fig1.tight_layout()


#%% FIGURE 2: t vs x   t vs y
# plt.rcParams['font.size'] = 10
# plt.rcParams["figure.figsize"] = (6,3)
# fig2, axes = plt.subplots(nrows=1, ncols=2)
    
# C = 0   
# axes[C].set_xlabel('t',color= 'black',fontsize = 12)
# axes[C].set_ylabel('x',color = 'black',fontsize = 12)
# axes[C].xaxis.grid()
# axes[C].yaxis.grid()
# axes[C].plot(t, xS,lw = 2)

# C = 1   
# axes[C].set_xlabel('t',color= 'black',fontsize = 12)
# axes[C].set_ylabel('y',color = 'black',fontsize = 12)
# axes[C].xaxis.grid()
# axes[C].yaxis.grid()
# axes[C].plot(t, yS,lw = 2)
# fig2.suptitle('$\Delta$t = %0.1f' %t2,color = 'black',fontsize = 14)
# fig2.tight_layout()


# #%% FIGURE 4: Phase Portrait  quiver plot
# plt.rcParams['font.size'] = 10
# plt.rcParams["figure.figsize"] = (5,4)
# fig4, axes = plt.subplots(nrows=1, ncols=1)
  
# axes.set_xlabel('x',color= 'black',fontsize = 12)
# axes.set_ylabel('y',color = 'black',fontsize = 12)
# axes.set_title('$\Delta$t = %0.1f' %t2,color = 'black',fontsize = 14)
# axes.xaxis.grid()
# axes.yaxis.grid()
# axes.set_xlim([L1,L2])
# axes.set_ylim([L1,L2])

# axes.quiver(xx,yy,xxDot,yyDot)

# axes.plot(xS,yS,'g',lw = 2)
# axes.plot(x0,y0,'go', ms = 6)

# axes.plot(xC,yC,'ro', ms = 7)

# axes.plot(xN0,yN0,'r', lw = 1)
# axes.plot(xN1,yN1,'b', lw = 1)
# axes.plot(xN1,-yN1,'b', lw = 1)
# axes.set_box_aspect(1)

# fig4.tight_layout()

# #%% FIGURE 5: Phase Portrait  streamine
# plt.rcParams['font.size'] = 10
# plt.rcParams["figure.figsize"] = (5,4)
# fig5, axes = plt.subplots(nrows=1, ncols=1)
  
# axes.set_xlabel('x',color= 'black',fontsize = 12)
# axes.set_ylabel('y',color = 'black',fontsize = 12)
# axes.set_title('$\Delta$t = %0.1f' %t2,color = 'black',fontsize = 14)
# axes.xaxis.grid()
# axes.yaxis.grid()
# axes.set_xlim([L1,L2])
# axes.set_ylim([L1,L2])
# axes.streamplot(xx,yy,xxDot,yyDot)



# axes.plot(xC,yC,'ro', ms = 7)

# axes.plot(xN0,yN0,'r', lw = 1)
# axes.plot(xN1,yN1,'b', lw = 1)
# axes.plot(xN1,-yN1,'b', lw = 1)

# axes.plot(xS,yS,'g',lw = 3)
# axes.plot(x0,y0,'go', ms = 6)

# axes.set_box_aspect(1)

# fig5.tight_layout()

# #%%   FIGURE 3: slope plot dy/dx
# plt.rcParams['font.size'] = 12
# plt.rcParams["figure.figsize"] = (5,4)
# fig3, ax = plt.subplots(nrows=1, ncols=1)
# XP = linspace(-2,2,255); YP = XP
# XX,YY = np.meshgrid(XP,XP)
# ZZ = (XX**2+YY**2-1)/(XX+1e-16)
# ZZ[(ZZ)>20]=20;ZZ[(ZZ)<-20] = -20  
# ax.set_xlabel('x',color= 'black',fontsize = 12)
# ax.set_ylabel('y',color = 'black',fontsize = 12)
# ax.set_title('slope function  dy/dx',color = 'black',fontsize = 14)
# ax.set_xticks([-2,-1,0,1,2])
# ax.set_yticks([-2,-1,0,1,2])
# cf = ax.pcolor(XX,YY,abs(ZZ)**0.3, cmap='jet') 
# ax.set_box_aspect(1)
# fig3.colorbar(cf, ax=ax)
# fig3.tight_layout()


# #%%
# tExe = time.time() - tStart
# print('  ')
# print('Execution time')
# print(tExe)

# '''
# #%%
'''
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png')
fig5.savefig('a5.png')


'''
