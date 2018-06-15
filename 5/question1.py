# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 16:16:39 2017

@author: Prabha
"""
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    a = (x*np.sin(x)-(np.sin(x)*np.sin(x)))/((1-np.cos(x))*(1-np.cos(x)))
    f = a + 2.182
    return f

diff = []

for i in np.arange(4,5,0.001):
    if f(i) >= 0:
        diff.append(i)
        
print (diff[-1])