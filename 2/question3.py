# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 13:40:19 2017

@author: Prabha
"""

#3a
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize as op

def f(x):
    y = ((np.exp(-x))*(x + 1)) - 0.5
    return y

x1 = 0.0
x2 = 2.0

print(f(x1), f(x2))

root = op.bisect(f,x1,x2)

print(root)


#3c

h = 5.4
I0 = 1.38e8

def L(x):
    y = 2*np.pi*I0*h*(h-(np.exp(-x/h)*(x+h)))
    return y

r = np.arange(0,100,0.01)

plt.figure()
plt.plot(r,L(r))
plt.title('Total Luminosity as a Function of Distance')
plt.xlabel('Distance r (kpc)')
plt.ylabel('Luminosity L (L$_{\odot}$)')
plt.savefig('3c.png', dpi = 1000)

