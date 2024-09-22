# -*- coding: utf-8 -*-
"""
qmHL.py
"""

# LIBRARIES  ==============================================================
import numpy
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, radians
import scipy.special 
import scipy.special as sps
from scipy import *
import scipy.linalg as la
import math
import cmath
import os
import matplotlib.pyplot as plt
from pylab import *
from math import *
from numpy import linspace,matrix,array
from numpy import *
from scipy.linalg import *
from scipy.optimize import fsolve
from scipy.integrate import odeint, solve_ivp, simps
from scipy.sparse.linalg import eigsh, eigs #Solves the Eigenvalue problem
from scipy.sparse import diags #Allows us to construct our matrices
from matplotlib.animation import FuncAnimation, PillowWriter 
import time
from scipy.special import sph_harm

tStart = time.time()


N = 599
t = linspace(0,pi,N)

c32 = -sqrt(21/(64*pi))*(sqrt(2*pi))

c20 = sqrt(5/(16*pi))*(sqrt(2*pi))

T31 = c32*(4*cos(t)**2*sin(t) - sin(t)**3)

T20 = c20*(3*cos(t)**2 - 1)

THETA32 = sph_harm(0, 0, 0, t).real

plt.plot(t,T31)
plt.plot(t,T20)


# Check normalization
fn = T31**2*sin(t)
A32 = simps(fn,t)
fn = T20**2*sin(t)
A20 = simps(fn,t)

print(A32,A20)
