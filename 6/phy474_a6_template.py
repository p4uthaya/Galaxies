import numpy as np
import matplotlib.pyplot as plt
import scipy
from astropy.io import ascii,fits
from astropy.table import Table
import astropy.coordinates as coordinates
import astropy.units as u
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator,DictFormatter)

def readTable_v1():
	# Here's the simplest way to read the data tables
	data=fits.getdata('PHY474_a6_cfaslice_dr14_mbalogh.fit')
	print (data.columns)
	# You can access the data now via the column name.  e.g. data['ra'] is the Right Ascension (in degrees)

def readTable_v2():
	# If you want to be slick about it, you can use astropy.coordinates and astropy.units to gracefully handle unit conversions (e.g. between radians and degrees).  
	t=Table.read('PHY474_a6_cfaslice_dr14_mbalogh.fit',format='fits')
	# Assign units to the RA and DEC columns.
	t['ra'].unit=u.deg
	t['dec'].unit=u.deg
	# Create a new field for the array t containing a coordinate object
	coordarr=coordinates.SkyCoord(t['ra'].quantity,t['dec'].quantity,frame='fk5')
	t['coord']=coordarr
	# you can now access ra and dec in different units, like this:
	print (t['coord'].dec.degree[1:10])
	print (t['coord'].dec.radian[1:10])
	#
def setup_axes3(fig, rect,ra0,ra1,z0,z1):
    """
    Sometimes, things like axis_direction need to be adjusted.
    """

    # rotate a bit for better orientation
    tr_rotate = Affine2D().translate(-100, 0)

    # scale degree to radians
    tr_scale = Affine2D().scale(np.pi/180., 1.)

    tr = tr_rotate + tr_scale + PolarAxes.PolarTransform()

    grid_locator1 = angle_helper.LocatorHMS(4)
    tick_formatter1 = angle_helper.FormatterHMS()

    grid_locator2 = MaxNLocator(3)

    grid_helper = floating_axes.GridHelperCurveLinear(
        tr, extremes=(ra0, ra1, z0, z1),
        grid_locator1=grid_locator1,
        grid_locator2=grid_locator2,
        tick_formatter1=tick_formatter1,
        tick_formatter2=None)

    ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
    fig.add_subplot(ax1)

    # adjust axis
    ax1.axis["left"].set_axis_direction("bottom")
    ax1.axis["right"].set_axis_direction("top")

    ax1.axis["bottom"].set_visible(False)
    ax1.axis["top"].set_axis_direction("bottom")
    ax1.axis["top"].toggle(ticklabels=True, label=True)
    ax1.axis["top"].major_ticklabels.set_axis_direction("top")
    ax1.axis["top"].label.set_axis_direction("top")

    ax1.axis["left"].label.set_text(r"Redshift")
    ax1.axis["top"].label.set_text(r"Right Ascension")

    # create a parasite axes whose transData in RA, cz
    aux_ax = ax1.get_aux_axes(tr)

    aux_ax.patch = ax1.patch  # for aux_ax to have a clip path as in ax
    ax1.patch.zorder = 0.9  # but this has a side effect that the patch is
    # drawn twice, and possibly over some other
    # artists. So, we decrease the zorder a bit to
    # prevent this.

    return ax1, aux_ax
def polar_example():
	fig=plt.figure(figsize=(10,10))
	# set RA limits in degrees
	ra0, ra1 = 130,250
	z0, z1 = 0, 0.1
	N=1000
	ra=np.random.uniform(size=N)*(ra1-ra0)+ra0
	z=np.random.uniform(size=N)*(z1-z0)+z0
	fig=plt.figure(figsize=(10,10))
	ax, aux_ax = setup_axes3(fig, 111,ra0,ra1,z0,z1)
	aux_ax.scatter(ra,z,s=1,color='SteelBlue',edgecolors="none",alpha=0.7)
	# The numpy 'where' function is very useful for plotting subsets of the data:
	colour=np.random.uniform(size=N)
	select=np.where(colour>0.5)
	aux_ax.scatter(ra[select],z[select],s=20,color='Red')

#1a

data=fits.getdata('PHY474_a6_cfaslice_dr14_mbalogh.fit')
print(data.columns)


u = np.array(data['u'])
r = np.array(data['r'])

ur = u-r
print (np.amax(ur))
print(np.amin(ur))
print (len(ur))
plt.figure()
plt.hist(ur,500)
plt.title('u-r Colour Magnitude Distribution')
plt.xlabel('u-r')
plt.ylabel('Number of Galaxies')
#plt.savefig('1a.png', dpi = 1000)

# 1b

ra0, ra1 = 120,262.4
z0, z1 = 0, 0.1

ra = np.array(data['ra'])
dec = np.array(data['dec'])
z = np.array(data['redshift'])

print(np.amax(z))

fig=plt.figure(figsize=(10,10))
ax, aux_ax = setup_axes3(fig, 111,ra0,ra1,z0,z1)
colour = ur
selectb=np.where(colour<2.5)
aux_ax.scatter(ra[selectb],z[selectb],s=1,color='SteelBlue',edgecolors="none",alpha=0.7)
selectr=np.where(colour>2.5)
aux_ax.scatter(ra[selectr],z[selectr],s=1,color='Red',edgecolors="none",alpha=0.7)


abdata=fits.getdata('abell_phy474.fits')
print(abdata.columns)

abra = np.array(abdata['RA'])
abz = np.array(abdata['REDSHIFT'])
abrch = abdata['RICHNESS']

sc =aux_ax.scatter(abra,abz,s=20,c=abrch,cmap='plasma',edgecolors="none",alpha=0.9)
cbar=plt.colorbar(sc,orientation='horizontal')
cbar.set_label('Richness')
#plt.savefig('1b1.png', dpi = 1000)


#1c

abname = abdata['NAME']
abdec = abdata['DEC']

CBgal= ['ACO 2056','ACO 2061','ACO 2065','ACO 2067','ACO 2079','ACO 2089','ACO 2092']


index = []
for i in range(len(abname)):
    for j in range(len(CBgal)):
        if abname[i] == CBgal[j]:
            index.append(i)

print(index)

CBra = []
CBdec = []
for i in index:
    CBra.append(abra[i])
    CBdec.append(abdec[i])

print(sum(CBra)/7)
print(sum(CBdec)/7)



rac = ra - 231.4
decc = dec - 29.4
CBrac = np.array(CBra)-231.4
CBdecc = np.array(CBdec)-29.4


plt.figure()

col = ur
selectb = np.where(col<2.5)
plt.scatter(rac[selectb],decc[selectb],s=1,color='SteelBlue')
selectr = np.where(col>2.5)
plt.scatter(rac[selectr],decc[selectr],s=1,color='Red')
plt.scatter(CBrac,CBdecc,alpha = 0.7,color='limegreen',edgecolors="none", s = 10)
axes = plt.gca()
axes.set_xlim([-3,3])
axes.set_ylim([-3,3])
plt.title('Corona Borealis Super Cluster')
plt.xlabel('Right Ascension from Centre ($\degree$)')
plt.ylabel('Declination from Centre ($\degree$)')
#plt.savefig('1b2.png', dpi = 1000)

c = 3e5
H = 70
abz = abdata['REDSHIFT']
CBz = []
for i in index:
    CBz.append(abz[i])

CBz = np.array(CBz)
CBra = np.array(CBra)
CBdec = np.array(CBdec)

print(CBz)
print(CBra)
print(CBdec)
print (z)
print (ra)
print (dec)


raoff = []
decoff = []
zoff = []
indexoff = []

for i in range(len(ra)):
    if 228.4 < ra[i] < 234.4:
        if 26.4 < dec[i] < 32.4:
            raoff.append(ra[i])
            decoff.append(dec[i])
            zoff.append(z[i])
            indexoff.append(i)

raoff = np.array(raoff)
decoff = np.array(decoff)
zoff = np.array(zoff)

for i in range(len(CBz)):
    vel = -c*(CBz[i]-zoff)
    delra = CBra[i] - raoff
    deldec = CBdec[i] - decoff
    delalpha = delra*np.cos(decoff)
    theta = np.sqrt(deldec*deldec + delalpha*delalpha)
    theta = (2*np.pi/360)*theta
    dc = (c/H)*(CBz[i])
    dis = dc*theta
    plt.figure()
    plt.scatter(dis,vel,s=1)
    plt.title('%s'%CBgal[i])
    axes = plt.gca()
    axes.set_xlim([0,10])
    plt.xlabel('Distance from Centre (Mpc)')
    plt.ylabel('Offset Velocity (km/s)')
    plt.savefig('1b3 %s .png'%CBgal[i], dpi = 100, bbox_inches='tight' )
    meanv = sum(vel)/len(vel)
    sddis = np.std(dis)
    print ('deviation for %s'%CBgal[i], sddis)















