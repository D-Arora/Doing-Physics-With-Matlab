# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 18:55:55 2024

@author: Owner
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 13:03:24 2020

@author: Mohammad Asif Zaman

Code for testng differentiation matrix operators

"""


#from __future__ import print_function    

import time
import math
import numpy as np
import pylab as py
import scipy.sparse as sp                 # import sparse matrix library

py.rcParams.update({'font.size': 14})

# import the file where the differentiation matrix operators are defined
#from diff_matrices import Diff_mat_1D, Diff_mat_2D   

def my_contourf(x,y,F,ttl):
    py.contourf(x,y,F,41,cmap = 'inferno')
    py.colorbar()
    py.xlabel('x'); py.ylabel('y'); py.title(ttl)
    return 0
    
def Diff_mat_1D(Nx):
    
    # First derivative
    D_1d = sp.diags([-1, 1], [-1, 1], shape = (Nx,Nx)) # A division by (2*dx) is required later.
    D_1d = sp.lil_matrix(D_1d)
    D_1d[0,[0,1,2]] = [-3, 4, -1]               # this is 2nd order forward difference (2*dx division is required)
    D_1d[Nx-1,[Nx-3, Nx-2, Nx-1]] = [1, -4, 3]  # this is 2nd order backward difference (2*dx division is required)
    
    # Second derivative
    D2_1d =  sp.diags([1, -2, 1], [-1,0,1], shape = (Nx, Nx)) # division by dx^2 required
    D2_1d = sp.lil_matrix(D2_1d)                  
    D2_1d[0,[0,1,2,3]] = [2, -5, 4, -1]                    # this is 2nd order forward difference. division by dx^2 required. 
    D2_1d[Nx-1,[Nx-4, Nx-3, Nx-2, Nx-1]] = [-1, 4, -5, 2]  # this is 2nd order backward difference. division by dx^2 required.
    
    return D_1d, D2_1d




def Diff_mat_2D(Nx,Ny):
    # 1D differentiation matrices
    Dx_1d, D2x_1d = Diff_mat_1D(Nx)
    Dy_1d, D2y_1d = Diff_mat_1D(Ny)


    # Sparse identity matrices
    Ix = sp.eye(Nx)
    Iy = sp.eye(Ny)


    
    # 2D matrix operators from 1D operators using kronecker product
    # First partial derivatives
    Dx_2d = sp.kron(Iy,Dx_1d)
    Dy_2d = sp.kron(Dy_1d,Ix)
    
    # Second partial derivatives
    D2x_2d = sp.kron(Iy,D2x_1d)
    D2y_2d = sp.kron(D2y_1d,Ix)
    
   
    
    # Return compressed Sparse Row format of the sparse matrices
    return Dx_2d.tocsr(), Dy_2d.tocsr(), D2x_2d.tocsr(), D2y_2d.tocsr()




# Define independent variables
Nx = 80                         # No. of grid points along x direction
Ny = 70                         # No. of grid points along y direction
x = np.linspace(-3,3,Nx)        # x variables in 1D
y = np.linspace(-3,3,Ny)        # y variable in 1D

dx = x[1] - x[0]                # grid spacing along x direction
dy = y[1] - y[0]                # grid spacing along y direction

X,Y = np.meshgrid(x,y)          # 2D meshgrid


# Test funcitons
F = np.exp(-X**2-3*Y**2)        # 2D Gaussian test function 
Fx = np.multiply(X,F)*(-2)      # Analytical value of dF/dx
Fy = np.multiply(Y,F)*(-6)      # Analytical value of dF/dy

F2x = np.multiply(2*X**2-1,F)*(2)      # Analytical value of d2F/dx2
F2y = np.multiply(6*Y**2-1,F)*(6)      # Analytical value of d2F/dy2


G = np.exp(-x**2)
Gx = np.multiply(x,G)*(-2)
G2x = np.multiply(2*x**2-1,G)*(2)


Dx_1d, D2x_1d = Diff_mat_1D(Nx)                     # Calling 1D matrix operators from function
Dx_2d, Dy_2d, D2x_2d, D2y_2d = Diff_mat_2D(Nx,Ny)   # Calling 2D matrix operators from funciton

# Differentiation using matrix operators. The division by 2dx and 2dy are needed as they are not defined within the funcitons

dGx =  Dx_1d*G/(2*dx)     # dF/dx from matrix operator
d2Gx =  D2x_1d*G/(dx**2)     # dF/dx from matrix operator

dFx = (Dx_2d*F.ravel()).reshape([Ny,Nx])/(2*dx)     # dF/dx from matrix operator
dFy = (Dy_2d*F.ravel()).reshape([Ny,Nx])/(2*dy)     # dF/dy from matrix operator

d2Fx = (D2x_2d*F.ravel()).reshape([Ny,Nx])/(dx**2)     # d2F/dx2 from matrix operator
d2Fy = (D2y_2d*F.ravel()).reshape([Ny,Nx])/(dy**2)     # d2F/dy2 from matrix operator

# Plotting 1D differntiations
py.close('all')
#py.contourf(x,y,F,41,cmap = 'plasma')

py.figure(figsize = (14,4.8))
py.subplot(1,2,1)
py.plot(x,Gx,'o',label = 'Analytical')
py.plot(x,dGx,label = 'Matrix')
py.legend()
py.xlabel('x'); py.ylabel(r'$dG/dx$'); py.title('First derivative')


py.subplot(1,2,2)
py.plot(x,G2x,'o',label = 'Analytical')
py.plot(x,d2Gx,label = 'Matrix')
py.legend()
py.xlabel('x'); py.ylabel(r'$d^2G/dx^2$'); py.title('Second derivative') 



# Plotting 2D differentiations

# First partial derivatives
py.figure(figsize = (14,4.8))
py.subplot(1,2,1); my_contourf(x,y,Fx,r'Analytical $dF/dx$')
py.subplot(1,2,2); my_contourf(x,y,dFx,r'Matrix operator $dF/dx$')


py.figure(figsize = (14,4.8))
py.subplot(1,2,1); my_contourf(x,y,Fy,r'Analytical $dF/dy$')
py.subplot(1,2,2); my_contourf(x,y,dFy,r'Matrix operator $dF/dy$')


# Second partial derivatives
py.figure(figsize = (14,4.8))
py.subplot(1,2,1); my_contourf(x,y,F2x,r'Analytical $d^2F/dx^2$')
py.subplot(1,2,2); my_contourf(x,y,d2Fx,r'Matrix operator $d^2F/dx^2$')


py.figure(figsize = (14,4.8))
py.subplot(1,2,1); my_contourf(x,y,F2y,r'Analytical $d^2F/dy^2$')
py.subplot(1,2,2); my_contourf(x,y,d2Fy,r'Matrix operator $d^2F/dy^2$')