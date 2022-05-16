"""
This is similar to BackgroundNRL.py but uses background computed using
the CIMP Background class.
"""

from matplotlib.colors import Colormap
from CIMP import Enhance

import math
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import sunpy.visualization.colormaps as cm

from sunpy.net import attrs as a
from sunpy.net import Fido
from sunpy.io import fits

import sunkit_image.utils.utils

from scipy import interpolate

#======================================================================
def radial_cut(r,a,N=100):
    """
    extract a curve vs theta at a selected radius
    radius is given in terms of the size of the image, assuming
    the range of the image goes from -1 to 1
    """

    ny, nx = data.shape[:2]

    x = np.arange(nx)
    y = np.arange(ny)

    theta = np.linspace(0.,2*np.pi,N)

    ainterp = interpolate.interp2d(x, y, a, kind='linear')    

    acut = np.zeros(N)

    for idx in np.arange(N):
        xp = (1.0 + r1 * np.cos(theta[idx])) * nx/2
        yp = (1.0 + r1 * np.sin(theta[idx])) * ny/2
    
        ap = ainterp(xp,yp)
    
        acut[idx] = ap[0]

    theta_deg = theta * 180.0 / np.pi

    return acut, theta_deg

#======================================================================
pcase = 1

if pcase == 1:

    dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04'
    file = dir+'/15/32296650.fts'
    bgfile = dir+'/'+'background.fts'

else:
   print("Pick a valid test case")
   exit()

#======================================================================
# basic image to work with

data, header = fits.read(file)[0]
amap = sunpy.map.Map(data,header)
print(f"data range: {amap.min()}, {amap.max()}")

# Background file computed by CIMP
bkg, bheader = sunpy.io.fits.read(bgfile)[0]
bmap = sunpy.map.Map(bkg,header)
print(f"bkg range: {bmap.min()}, {bmap.max()}")

#======================================================================
# Now plot
cmap = plt.get_cmap('stereocor2')

fig = plt.figure(figsize=[16,12])

ax = fig.add_subplot(2,3,1,projection=amap)
#amap.plot(clip_interval=[10,90]*u.percent,cmap=cmap)
amap.plot(cmap=cmap)

ax = fig.add_subplot(2,3,4,projection=bmap)
#bmap.plot(clip_interval=[10,90]*u.percent,cmap=cmap)
bmap.plot(cmap=cmap)

#ax = fig.add_subplot(2,3,3,projection=dmap)
#dmap.plot(clip_interval=[10,90]*u.percent)
#dmap.plot(clip_interval=[10,90]*u.percent)

#======================================================================
# polar intensity plots at a particular radius

r1 = 0.2
acut, theta = radial_cut(r1,data)
bcut, theta = radial_cut(r1,bkg)

ax = fig.add_subplot(2,3,2)
plt.plot(theta,acut,'black')
plt.plot(theta,bcut,'blue')
plt.title("r = 0.2")

r1 = 0.8
acut, theta = radial_cut(r1,data)
bcut, theta = radial_cut(r1,bkg)

ax = fig.add_subplot(2,3,5)
plt.plot(theta,acut,'black')
plt.plot(theta,bcut,'blue')
plt.title("r = 0.8")


#======================================================================
# difference

pa = data - bkg
pmap = sunpy.map.Map(pa,header)

ax = fig.add_subplot(2,3,3,projection=pmap)
#pmap.plot(clip_interval=[20,90]*u.percent,cmap=cmap)
pmap.plot(cmap=cmap)

#======================================================================
# ratio

rat = np.where(bkg <= 0.0, 0.0, data/bkg)
rmap = sunpy.map.Map(rat,header)

#emap = Enhance.fnrgf(rmap, instrument, detector)

ax = fig.add_subplot(2,3,6,projection=pmap)
#pmap.plot(clip_interval=[20,90]*u.percent,cmap=cmap)
pmap.plot(cmap=cmap)

plt.show()