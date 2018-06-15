# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 12:49:54 2017

@author: Prabha
"""

# question 3
import numpy as np
import matplotlib.pyplot as plt

# b)
Zsun = 0.0134
p = (Zsun*(1-0.15))/(np.log(50/13))

print (p)
print (p/Zsun)

Mstarfin=50*(1-np.exp((0.15*Zsun-0.25*Zsun)/p))

print(Mstarfin)

# c)

p = Zsun

z_list1 = []
z_list2 = []
z_list3 = []
frac_list =[]
frac_list2 = []

for i in np.arange(0.0001,1.0001,0.0001):
    z1 = (p/-0.5)*(1 - i**(0.5/(-0.5-1)))
    z2 = p*(np.log(1/i))
    z3 = (p/0.5)*(1 - i**(-0.5/(0.5-1)))

    z_list1.append(z1)
    z_list2.append(z2)
    z_list3.append(z3)
    frac_list.append(i)

plt.figure()
plt.plot(frac_list, z_list1)
plt.title('Metallicity when v = -0.5')
plt.xlabel('$M_g$(t)/$M_g$(0)')
plt.ylabel('Metallicity Z(t)')
plt.savefig('3c1.png', dpi = 1000)

plt.figure()
plt.plot(frac_list, z_list2)
plt.title('Metallicity when v = 0')
plt.xlabel('$M_g$(t)/$M_g$(0)')
plt.ylabel('Metallicity Z(t)')
plt.savefig('3c2.png', dpi = 1000)

plt.figure()
plt.plot(frac_list, z_list3)
plt.title('Metallicity when v = 0.5')
plt.xlabel('$M_g$(t)/$M_g$(0)')
plt.ylabel('Metallicity Z(t)')
plt.savefig('3c3.png', dpi = 1000)