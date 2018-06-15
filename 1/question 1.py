# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import random as rand
import matplotlib.pyplot as plt

# Question 1

# 1a) 
i=2000
rmin= 70
rmax= 130
pi=np.pi

dVmin = round (4*pi*(rmin)**2)
dVmax = round (4*pi*(rmax)**2)

# print(dVmin,dVmax) - checking the range

r_list = []

for i in range(2000):
    V = rand.randint(dVmin,dVmax)
    r = np.sqrt(V/(4*pi))
    r_list.append(r)
    
# print(len(r_list)) - checking the number of values
# print(r_list[0]) - checkinh the order of values

Mvsun = 4.83
sig = 0.3

mag_list = np.random.normal(Mvsun,sig,2000)

plt.plot(r_list,mag_list,'.')
plt.title('Absolute magnitudes of G-type Stars as a Function of Distance')
plt.xlabel('Distance r (pc)')
plt.ylabel('Absolute Magnitude $M_v$')


# 1b) 
Mvis_list = []
rvis_list = []
apmagvis_list= []
apmagall_list= []

for i in range(2000):
    apmag = mag_list[i] + 5*np.log10(r_list[i]/10)
    apmagall_list.append(apmag)
    if apmag < 10:
        Mvis_list.append(mag_list[i])
        rvis_list.append(r_list[i])
        apmagvis_list.append(apmag)
    

#print (len(Mvis_list))
#print (len(rvis_list))
#print(len(apmagvis_list)) 
#print(len(apmagall_list))

meanmag = sum(Mvis_list)/len(Mvis_list)
print ("mean absloute magnitude of the sample of detected stars: %s" % meanmag)

plt.figure(1)
vis = plt.plot(rvis_list,Mvis_list,'.')
plt.savefig('1ab.png', dpi = 1000)

# 1c)
avgr = sum(rvis_list)/len(rvis_list)
print("average distance of the sample of detected stars: %s" % avgr)

L_list = []

msun = -26.73
rsun = 4.848e-6

for i in range(2000):
    L = (10**((msun-apmagall_list[i])/2.5))*((r_list[i])**2/rsun**2)
    L_list.append(L)

avgL = sum(L_list)/len(L_list)
print("average L for all stars is %s solar luminosities" % avgL)

rnew_list = []

samp = len(apmagvis_list)
for i in range(samp):
    rnew = np.sqrt((avgL)/(10**((msun-apmagvis_list[i])/2.5)))*rsun
    rnew_list.append(rnew)

avgrnew = sum(rnew_list)/len(rnew_list)
print("average distance of the sample of detected stars with Lavg : %s" % avgrnew)

#1d)

zmean = 0.67
zsig = 0.28

z_list = np.random.normal(zmean,zsig,2000)

for i in range(2000):
    if z_list[i] <= 0. :
        z_list[i] = 0.0000000001

Mnew_list = []
zvis_list = []
rzvis_list = []

for i in range(2000):
    delM = -0.87*np.log10(z_list[i])
    Mnew = mag_list[i]+delM
    Mnew_list.append(Mnew)
    mnew = Mnew_list[i] + 5*np.log10(r_list[i]/10)
    if mnew < 10:
        zvis_list.append(z_list[i])
        rzvis_list.append(r_list[i])

plt.figure(2)
plt.plot(rzvis_list,zvis_list,'.')
plt.show

bin1 = []
bin2 = []
bin3 = []
bin4 = []
bin5 = []
bin6 = []

for i in range(len(rzvis_list)):
    if  70 <= rzvis_list[i] < 80:
        bin1.append(zvis_list[i])
    elif 80 <= rzvis_list[i] < 90:
        bin2.append(zvis_list[i])
    elif 90 <= rzvis_list[i] < 100:
        bin3.append(zvis_list[i])
    elif 100 <= rzvis_list[i] < 110:
        bin4.append(zvis_list[i])
    elif 110 <= rzvis_list[i] < 120:
        bin5.append(zvis_list[i])
    elif 120 <= rzvis_list[i] <= 130:
        bin6.append(zvis_list[i])

avgmet = []
binavg1 = sum(bin1)/len(bin1)
avgmet.append(binavg1)

binavg2 = sum(bin2)/len(bin2)
avgmet.append(binavg2)

binavg3 = sum(bin3)/len(bin3)
avgmet.append(binavg3)

binavg4 = sum(bin4)/len(bin4)
avgmet.append(binavg4)

binavg5 = sum(bin5)/len(bin5)
avgmet.append(binavg5)

binavg6 = sum(bin6)/len(bin6)
avgmet.append(binavg6)

print(avgmet)

plt.plot([75,85,95,105,115,125],avgmet)
plt.title('Metallicities of G-type Stars as a Function of Distance')
plt.xlabel('Distance r (pc)')
plt.ylabel('Metallicity $Z/Z_\odot$')
plt.savefig('1d.png', dpi = 1000)





    


