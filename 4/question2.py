# -*- coding: utf-8 -*-
"""
Created on Sun Nov  5 20:11:00 2017

@author: Prabha
"""
import numpy as np
import matplotlib.pyplot as plt

Msun = 1.99e30
Mbh = 4e6*Msun
ro = 4e6*Msun #pc^-3
a = 1         #pc
G = 6.67e-11*(1/3.086e16)*(1/3.086e16)*(1/3.086e16)

# luminous stars only

def vd(r):
    x = r/a
    vdsq = (2/9)*np.pi*G*ro*a*a*(1/np.sqrt(1+(x*x)))
    vd = np.sqrt(vdsq)*3.086e13
    return vd

r = np.arange(0.1,10,0.1)

plt.figure()
plt.plot(r,vd(r))
plt.title('Velocity Dispersion - Star Cluster')
plt.xlabel('Distance r (pc)')
plt.ylabel('Velocity Dispersion  $\sigma$$_{r}$ (km/s)')
plt.savefig('2c1.png', dpi = 1000)

# stars + bh

def vdb(r):
    x = r/a
    b = ((G*Mbh)/a)*((np.sqrt(1+x*x))**5)*((-8/3)+((8*x**4 + 12*x**2 + 3)/(3*x*(np.sqrt(1+x*x))**3)))
    vdsq = (2/9)*np.pi*G*ro*a*a*(1/np.sqrt(1+(x*x))) + b
    vdb = np.sqrt(vdsq)*3.086e13
    return vdb

def sc(r):
    x = r/a
    vdsq = (2/9)*np.pi*G*ro*a*a*(1/np.sqrt(1+(x*x)))
    return np.sqrt(vdsq)*3.086e13

def bh(r):
    x = r/a
    b = ((G*Mbh)/a)*((np.sqrt(1+x*x))**5)*((-8/3)+((8*x**4 + 12*x**2 + 3)/(3*x*(np.sqrt(1+x*x))**3)))
    return np.sqrt(b)*3.086e13


plt.figure()
plt.plot(r,vdb(r))
plt.plot(r,sc(r),':')
plt.plot(r,bh(r),'--')
plt.title('Velocity Dispersion - Star Cluster with Black Hole')
plt.xlabel('Distance r (pc)')
plt.ylabel('Velocity Dispersion  $\sigma$$_{r}$ (km/s)')
plt.savefig('2c2.png', dpi = 1000)


diff = []

for i in np.arange(0.001,1,0.001):
    if bh(i) - sc(i) >= 0:
        diff.append(i)
        
print (diff[-1])
    