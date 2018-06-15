# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 23:30:16 2017

@author: Prabha
"""
# Question 2

# a)

import numpy as np
import matplotlib.pyplot as plt

msun = 1.989e30


frac_list = []
tau_list = []

for i in np.arange(2,10.01,0.01):
    NT = (1-np.exp(-10/i))
    NV = (np.exp(-8.9/i)-np.exp(-10/i))
    frac = NV/NT
    frac_list.append(frac)
    tau_list.append(i)

plt.plot(tau_list, frac_list)
plt.title('Fraction of Total 2$M_\odot$ Stars Visible on Main Sequence Today')
plt.xlabel('Timescale $\u03C4$ (Gyr)')
plt.ylabel('Fraction of Visible Stars')
plt.savefig('2a.png', dpi = 1000)
# b)

m1 = 2*msun
m2 = 0.8*msun

n1 = 1
n2 = 100

mf1 = n1/m1
mf2 = n2/m2

