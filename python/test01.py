# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 14:48:30 2022

@author: Owner
"""


# Reference: python_221017.pdf

import numpy

import numpy as np

# Special functions
# https://docs.scipy.org/doc/scipy/reference/special.html



import scipy.special 

import scipy.special as sps

import math

import os

import matplotlib.pyplot as plt


print('Hello world')
print('Ian first program')
a = 5
b = 7
c = a+b
print(c)

# age = input('your age   ')

# print(age)

x = 22
y = 5;
z1 = x+y
z2 = x-y
z3 = x*y
z4 = x/y
z5 = x//y
z6 = x%y
z7 = x**y

s = 'This is a string'
r = '  string functions'
s1 = len(s)
s2 = s.count('i')
s3 = s.lower()
s4 = s[5]
s5 = s[1:16:3]

print(s1)
print(s2)
print(s3)
print(s4)
print(s5)
print(s+r)


a = 72.455
b = 4536.66778

print(a+b)
print('The first number is {:2.0f} and the second number is {:3.2e}'.format(a, b))

joe_string='My GPA is {:5.2f} and I am {:5.3f} years old.'
print(joe_string.format(b,a))


z = round(b)
print(z)

x = 2.14
y = 2

b,c = divmod(x,y)


x = 0.5*numpy.pi # Get the value of pi
y = numpy.sin(x) # Find the sine of pi radians
z = numpy.sqrt(x)


a = np.pi

print(a)



a = 5.5
b = 6.2
c = np.pi/2
d = a * np.sin(b * c)
e = a * sps.j0(3.5 * np.pi/4)

print(e)

#%%   CELL
#os.system('cls')
# Arrays    Array Creation
# s=sum((a5)**3)
# In most respects, math with Numpy arrays behaves like the “dotted” operators
#   with Matlab matrices, where the operations are performed on corresponding
#   elements.

a = np.linspace(1,10,10)

s = sum((a-5)**3)
print(a,s)

print()

a = np.array([1,2,3,4,5,6,7,8,9,10])
print(a)

x = np.linspace(1,10,10)
y = np.linspace(3,23,10)
x1 = x**2
x2 = x+y
x3 = x + 3
x4 = x*y


a = x

print
print(x)
print(a)

a = np.linspace(3,23,10)

print()

print(x)
print(a)


#%% CELL     xy graphs


x = np.linspace(0,10,21)
y = x**2
y
#plt.figure(1)
#plt.plot(x,y,'or')

x = np.linspace(0,4*np.pi,210)
y1 = np.sin(3*x)
y2 = np.cos(2*x)


plt.figure(2)
plt.plot(x,y1,'k')
plt.plot(x,y2,'r')



plt.legend(['Sine of x','Cosine of x'])
plt.xlabel('x')
plt.ylabel('y')
plt.title('Two Trig Functions')
plt.grid(color='k', linestyle='-', linewidth = 0.2)

#plt.figure(3)
#plt.plot(y1,y2)
#plt.title('Something Fancy')

#%% CELL

s = 0
n = 0


#x = np.linspace(0,4*np.pi,5) 
n = np.linspace(1,10,10) 
print(n)
print('***')
for x in n:
    s += x
    print(s)
    
    
#for n in range(0,10):
 #   s += n
  #  print(s)
   
   
#%%
"""
Example showing shaded 3d plots. It is based on the [shading example](
http://matplotlib.org/examples/pylab_examples/shading_example.html).
The surface used is the Matlab `peaks()`.
"""

from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource

fig = plt.figure()
#ax = fig.gca(projection='3d')
ax = plt.subplots(subplot_kw={"projection": "3d"})
# Test data: Matlab `peaks()`
x, y = np.mgrid[-3:3:150j,-3:3:150j]
z =  3*(1 - x)**2 * np.exp(-x**2 - (y + 1)**2) \
   - 10*(x/5 - x**3 - y**5)*np.exp(-x**2 - y**2) \
   - 1./3*np.exp(-(x + 1)**2 - y**2) 

# create light source object.
ls = LightSource(azdeg=0, altdeg=65)
# shade data, creating an rgb array.
rgb = ls.shade(z, plt.cm.RdYlBu)
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, linewidth=0,
                       antialiased=False, facecolors=rgb)


plt.show()




