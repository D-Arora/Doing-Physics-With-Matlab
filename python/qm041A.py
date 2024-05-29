# -*- coding: utf-8 -*-
"""
Created on Fri May 10 19:59:52 2024

@author: Owner



"""

### IMPORTS ###
import numpy as np
from scipy.sparse import diags #Allows us to construct our matrices
from scipy.sparse.linalg import eigsh #Solves the Eigenvalue problem
import matplotlib.pyplot as plt

### CONSTANTS ###
N = 100
hbar = 1.054571800e-34         # Reduced Planck's constant
m = 9.10938356e-31            # Electron mass [kg]
xmin = -1e-9
xmax = 1e-9
x = np.linspace(xmin, xmax, N)
dx = (xmax - xmin)/(N + 1)
L = 1.5e-9
W = 1e-3
U = np.where((x < L/2) & (x > -L/2), 0, W)

k = 5
scale = 50

### PREPARE MATRICES ###
d2dx2 = diags([1, -2, 1], offsets=[-1, 0, 1], shape=(N, N))/(dx**2)
T = -(hbar**2)/(2 * m) * d2dx2
V = diags(U)
H = T + V

### CALCULATE EIGENSTATES ###
eigvals, eigvecs = eigsh(H, which="SM", k=k)

### PLOTTING ###
# plt.plot(x, U)
# for i in range(k):
#     plt.plot(x, np.full_like(x, eigvals[i]), "k--")
#     plt.plot(x, eigvals[i] + eigvecs[:, i] * scale, 'k')
#     plt.plot(x, eigvals[i] + eigvecs[:, i] * eigvecs[:, i].conj() * scale * scale, 'b')
# plt.show()

plt.plot(x,eigvecs[:,3])

