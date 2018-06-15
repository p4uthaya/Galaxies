# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 14:46:23 2017

@author: Prabha
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

# 3c
def f(x):
    t = ((np.log(x+1))/x) + (1/(x*(x+1))) - (1/x)
    tp = ((x-(x+1)*np.log(x+1))/(x*x*(x+1))) - ((2*x+1)/(x*x*(x+1)*(x+1))) + (1/(x*x))
    vp = t**(-0.5)*tp
    return vp
def v(x):
    t = ((np.log(x+1))/x) + (1/(x*(x+1))) - (1/x)
    v = t**(0.5)
    return v

r = optimize.bisect(f,1,3)
vmax = v(r)

print(r)
print(vmax)

# 3f

roc = 136                          # Msun/kpc^3
G = 6.674e-11                      # m^3 kg^-1 s^-2
Msun = 1.989e30                    # kg


def t(rv,c,dl):
    ro = dl*roc                    # density parameter Msun/kpc^3
    a = rv/c                       # length parameter kpc
    x = 2.1625815870629594
    r = x*a                        # r with velocity is max kpc
    Mn = 4*np.pi*ro*(a*a*a)        # Mass parameter Msun
    Vn = ((G*(3.404e-59)*Mn*Msun)/a)**(0.5) # Velocity parameter kpc/s
    
    M = Mn*((np.log(c+1)) + (1/(c+1)) - 1)
    dn = (3*dl*(a*a*a)*((np.log(c+1)) + (1/(c+1)) - 1))/(rv*rv*rv) 
    vmax = (Vn*(((np.log(x+1))/x) + (1/(x*(x+1))) - (1/x))**(0.5)) * (3.086e16)
    sig = (Vn*((1/6)*(Mn/M)*((c-np.log(c+1))/(c+1)))**(0.5)) * (3.086e16)
    
    return M, dn, vmax, 0.465*Vn*3.086e16, r, sig

print(t(96,9.0,34600))

print(t(206,7.14,20000))

print(t(2060,3.55,4050))





