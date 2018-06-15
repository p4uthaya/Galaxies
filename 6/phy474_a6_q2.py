# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 00:39:47 2017

@author: Prabha
"""
import numpy as np
import matplotlib.pyplot as plt

lsun = 3.826e26
lam = 3.74e-6
r = 58600*3.086e22
c = 3e8
jkc = 10e26

def fv(l):
    fv = (l*lsun*lam*jkc)/(4*np.pi*r*r*c)
    m = -2.5*np.log10(fv)+8.90
    return fv, m

print(fv(2e9))
    

