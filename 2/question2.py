import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy import optimize
from scipy import stats
from astropy import constants

def q2():
    data = ascii.read('phy474_ass2_SNdata.txt')
    x = data['log10cz']
    y = data['mmax']
    return x,y


# 2b

x = q2()[0]
y = q2()[1]

plt.figure()
plt.scatter(x,y)
plt.title('Apparent Magnitude at Maximum Light as a Function of log$_{10}$(cz)')
plt.xlabel('log$_{10}$(cz)')
plt.ylabel('Apparent Magnitude m$_{max}$')


#2c


def func(x,m,b):
    return(m*x + b)

fit = optimize.curve_fit(func,x,y)

m,b = fit[0]
error = fit[1]

merr = error[0][0]
berr = error[1][1]
print(m,b)
print(np.sqrt(merr), np.sqrt(berr))


#2d
def func2(x,b):
    return(5*x + b)

fixfit = optimize.curve_fit(func2,x,y)
b = fixfit[0][0]
berr = fixfit[1][0][0]
print(b, np.sqrt(berr))

yfit = func2(x,b)

print(yfit[1])
print(y[1])

plt.plot(x,yfit,'k')
plt.savefig('2bd.png', dpi = 1000)
siglist =[]
for i in range(0,12):
    ymu = y[i]-yfit[i]
    siglist.append(ymu**2)

print(siglist)

sigma = np.sqrt((1/12)*(sum(siglist)))
print(sigma)

    