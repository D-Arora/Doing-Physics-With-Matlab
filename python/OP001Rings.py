# -*- coding: utf-8 -*-

# op001.py

# Libraries
import numpy as np
from numpy import pi, sin, cos, linspace 
from numpy.linalg import eig
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

tStart = time.time()


def simpson1d(num):
    if num%2 == 0:
       num = num-1
    sc = np.ones(num)
    R = np.arange(1,num,2);   sc[R] = 4
    R = np.arange(2,num-1,2); sc[R] = 2
    return [num, sc]

# Ring structure
#   radius of ring r  [m]     no. of data points around a ring n
#   Greater the circumference of a ring -->  more grid points
#   Width of each ring  dr
#   Total no. grid points for Aperture  nQ
n1 = 100            # Aperture: grid points for inner ring
n2 = 500            # Aperture: grid points for outer ring
nR = 501            # Aperture: number of rings    must be ODD
a = 1          # radius of circular aperture  [m]
rMax = a
rMin = 1e-20
r = linspace(rMin, rMax, nR)
dr = r[1] - r[0]
m = (n2-n1) / (nR-1)
b = n2 - m * nR
ax =0; bx = 1;
Iring = np.zeros(nR)
for c in range(0,nR-1,1):
    num = 2*round(0.5*(m * c + b))+1
    num, sc  = simpson1d(num)
    h = (bx-ax)/(num-1)
    f = r[c]
    Iring[c] = f*sum(sc)

I = dr*(h/3)*sum(Iring)

print(I)


# >>>>> Input number of grid points: N must be odd
# N = 2999
# if N%2 == 0:
#     N = N-1
    
# # >>>>> Input function to be integrated
# #ax = 0; bx = 2; ay = 1; by = 5
# ax = 0; bx = 1; ay = 0; by = 1
# #ax = 0; bx = pi/2; ay = 0; by = pi/2
# #ax = -1; bx = 1; ay = -1; by = 1; a = 1
# x = linspace(ax,bx,N)
# y = linspace(ay,by,N,N)
# xx, yy = np.meshgrid(x,y)

# #f = xx**2*yy**3
# f = 6*xx**0*yy**0
# #f = cos(xx)*sin(yy)

# #f = np.real(a**2 - xx**2 - yy**2)
# #f[(xx**2 + yy**2) >= a**2] = 0
# #f = f**0.5

# f[yy > 1 - xx] = 0

# # Simpson [2D] coefficients
# sc = np.ones(N)
# R = np.arange(1,N,2);   sc[R] = 4;
# R = np.arange(2,N-1,2); sc[R] = 2
# scx, scy = np.meshgrid(sc,sc)
# sc2D = scx*scy

# # Calculate integral
# hx = (bx-ax)/(N-1); hy = (by-ay)/(N-1)
# h = hx * hy / 9
# integral = h*sum(sum(sc2D*f))

# print('\n Integral  =  ',integral)



# #%% GRAPHICS
# plt.rcParams['font.size'] = 10
# fig = plt.figure(figsize=(4,3))
# ax = plt.axes(projection='3d')
# ax.plot_surface(xx, yy, f, cmap='jet',
#       edgecolor='none', alpha=1,antialiased=True)
# ax.set_xlabel('x', fontsize=12)
# ax.set_ylabel('y', fontsize=12)
# ax.set_zlabel('f', fontsize=12)
# fig.tight_layout()

# # fig.savefig('a1.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)



#%%
from scipy import integrate 
gfg = lambda x: x**2 + x + 1
  
# using scipy.integrate.quad() method 
geek = integrate.quad(gfg, 1, 4) 
  
print(geek)

#%%
import numpy as np 
from scipy import integrate 
  
#x = np.arange(-4, 4,0.001)
x = np.linspace(-4,4,999) 
#y = np.sqrt(x) 
y = -4*x**4 + 20*x**3 - 40*x**2 - 320*x + 1664
# using scipy.integrate.simps() method 
gfg = integrate.simps(y, x) 
  
z = -4
I1 = -(4/5)*z**5 + (20/4)*z**4 - (40/3)*z**3 - (320/2)*z**2 + 1664*z
z = 4
I2 = -(4/5)*z**5 + (20/4)*z**4 - (40/3)*z**3 - (320/2)*z**2 + 1664*z

I = I2 - I1

print(gfg) 
print(I)


#%%
# def func(x, y):
#     if x == y == 0.0:
#         return 0.0
#     else:
#         return (1-np.sqrt(x**2+y**2))/(np.sqrt(x**2+y**2))*(x**2)
    
# x1 = -1; x2 = 1;
# y1 = -np.sqrt(1-x**2); y2 =-sqrt(1-x**2)    
# print (dblquad(func, x1, x2, y1, y2))

from scipy import pi, sqrt
from scipy.integrate import dblquad

def fxy(x,y):
    fxy = x**2 * y**3
    return fxy

x1 = 0; x2 = 2; y1 = 1; y2 = 5

y3 = linspace(1,5,57)
x3 = linspace(0,2,57)
I = (dblquad(fxy, 1, 5, 0, 2))
#I3 = (dblquad(fxy, y3, x3))
print(I)


func = lambda y, x: x**2 * y**3
print(dblquad(func, 0, 2, lambda x: 1, lambda x: 5))

z = fxy(1,2)
print(z)


#%%

def int_2D( x, y, xlo=0.0, xhi=10.0, ylo=0.0, yhi=10.0, Nx=100, Ny=100 ):

    # construct f(x,y) for given limits
    #-----------------------------------
    xlin = np.linspace(xlo, xhi, Nx)
    ylin = np.linspace(ylo, yhi, Ny)
    X,Y = np.meshgrid(xlin, ylin)

    f = 1.0 * np.exp(-((X - 8.)**2 + (Y - 8)**2)/(8.0))
    f += 0.5 * np.exp(-((X - 1)**2 + (Y - 9)**2)/(10.0))

    # construct 2-D integrand
    #-----------------------------------
    m = np.sqrt( (x - X)**2 + (y - Y)**2 )+0.001
    y = 1.0 / np.arcsinh( m ) * f

    # do a 1-D integral over every row
    #-----------------------------------
    I = np.zeros( Ny )
    for i in range(Ny):
        I[i] = np.trapz( y[i,:], xlin )

    # then an integral over the result
    #-----------------------------------    
    F = np.trapz( I, ylin )

    return F
x = linspace(0,10,100)
y = linspace(0,10,100)

F = int_2D(x,y)

print(F)
