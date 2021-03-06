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
# copied over from the Enhance class for convenience
def mask_annulus(im, rmin = 0.0, rmax = None, missingval = 0.0):
    """
    This sets the pixels inside rmin and/or outside rmax to the missing value (default 0)
    """

    nx = im.shape[0]
    ny = im.shape[1]

    nn = 0.5 * float(np.min((nx, ny)))

    rr1 = (rmin*nn)**2

    if rmax is None:
        rr2 = np.inf
    else:
        rr2 = (rmax*nn)**2

    x0 = 0.5*float(im.shape[0])
    y0 = 0.5*float(im.shape[1])

    for i in np.arange(0,im.shape[0]):
        for j in np.arange(0,im.shape[1]):
            r2 = (float(i)-x0)**2 + (float(j)-y0)**2
            if (r2 < rr1) or (r2 > rr2):
                im[i,j] = missingval

#======================================================================
pcase = 3

if pcase == 1:
    instrument = 'lasco'
    detector = 'c3'
    dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04'
    file = dir+'/15/32296650.fts'
    bgfile = dir+'/'+'background.fts'

elif pcase == 2:
    instrument = 'lasco'
    detector = 'c3'
    dir = '/home/mark.miesch/data/lasco_monthly/c3/2014_01'
    file = dir+'/17/33385479.fts'
    bgfile = dir+'/'+'background.fts'

elif pcase == 3:
    instrument = 'stereo'
    detector = 'cor2'
    dir = '/home/mark.miesch/sunpy/data/secchi_cor2/L1/2012/09'
    #file = dir+'/20/20120920_172400_14c2A.fts'
    file = dir+'/20/20120920_002400_14c2A.fts'
    bgfile = dir+'/'+'background.fts'

else:
   print("Pick a valid test case")
   exit()

#======================================================================
# basic image to work with
data, header = fits.read(file)[0]

# Background file computed by CIMP
bkg, bheader = sunpy.io.fits.read(bgfile)[0]

#======================================================================
# apply a mask fot data outside the FOV
pos_pix = np.ma.masked_less_equal(bkg, 0.0)
bgmin = pos_pix.min()
print(f"background missing val {bgmin}")

mask_annulus(data, rmin = 0.15, rmax = 1.0)
mask_annulus(bkg, rmin = 0.15, rmax = 1.0, missingval = bgmin)

#======================================================================
amap = sunpy.map.Map(data,header)
print(f"data range: {amap.min()}, {amap.max()}")

bmap = sunpy.map.Map(bkg,header)
print(f"bkg range: {bmap.min()}, {bmap.max()}")

#======================================================================
# Now plot
cmap = plt.get_cmap('stereocor2')

fig = plt.figure(figsize=[16,12])

ax = fig.add_subplot(2,3,1,projection=amap)
amap.plot(cmap=cmap,vmin=0,vmax=1.e-8)

ax = fig.add_subplot(2,3,4,projection=bmap)
bmap.plot(cmap=cmap,vmin=0,vmax=1.e-8)

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
pb = Enhance.powerlaw(pa)
pmap = sunpy.map.Map(pb,header)

#dmap = Enhance.nrgf(pmap, instrument, detector)

ax = fig.add_subplot(2,3,3,projection=pmap)

print(f"Difference minmax: {pmap.min()} {pmap.max()}")

#pmap.plot(cmap=cmap,vmin=0,vmax=1.0e3)
pmap.plot(cmap=cmap,vmin=0,vmax=1.0e-9)
#pmap.plot(cmap=cmap)
#dmap.plot(cmap=cmap,clip_interval=(1,99)*u.percent)

#======================================================================
# ratio

rat = np.where(bkg <= 0.0, 0.0, data/bkg)

rmap = sunpy.map.Map(rat,header)

#pa = Enhance.powerlaw(rat,n=1.0)
#rmap = sunpy.map.Map(pa,header)

#emap = Enhance.fnrgf(rmap, instrument, detector)
#emap = Enhance.nrgf(rmap, instrument, detector)

ax = fig.add_subplot(2,3,6,projection=rmap)

print(f"Ratio minmax: {rmap.min()} {rmap.max()}")

rmap.plot(cmap=cmap,vmin=1.0,vmax=3)

plt.show()