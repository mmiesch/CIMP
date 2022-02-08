"""
Experiment with different types of background images and see how they affect the generation of "pretties"
"""

from CIMP import Enhance

import math
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from sunpy.net import attrs as a
from sunpy.net import Fido
import sunpy.io
import sunpy.map

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
pcase = 2

thresh = None

if pcase == 1:

    instrument = 'lasco'
    detector = 'c3'
    file = '/home/mark.miesch/sunpy/data/lasco_c3/32305543.fts'
    bgfile = '/home/mark.miesch/data/sswdb/lasco/monthly/3m_clcl_120716.fts'
    thresh = 1.4

elif pcase == 2:
    instrument = 'lasco'
    detector = 'c3'
    file = '/home/mark.miesch/data/lasco_ssw/32305543.fts'
    bgfile = '/home/mark.miesch/data/lasco_ssw/32305543_bkg.fts'
    #thresh = 1.4

else:
    instrument = 'lasco'
    detector = 'c3'
    file = '/home/mark.miesch/data/lasco_ssw/32352249.fts'
    bgfile = '/home/mark.miesch/data/lasco_ssw/32352249_bkg.fts'
    thresh = 1.4


#======================================================================
# basic image to work with

data, header = sunpy.io.fits.read(file)[0]

if thresh is not None:

    # threshold to block out the occulter
    data = data.clip(min=thresh*np.min(data))
    vmask = thresh*np.min(data)
    a = np.where(data <= vmask, 0, data)
    data = a

amap = sunpy.map.Map(data,header)
print(f"data range: {amap.min()}, {amap.max()}")

# Background file provided by NRL

bkg, bheader = sunpy.io.fits.read(bgfile)[0]

if thresh is not None:
    bkg = bkg.clip(min=thresh*np.min(bkg))

bmap = sunpy.map.Map(bkg,header)
print(f"bkg range: {bmap.min()}, {bmap.max()}")

#======================================================================
# Now plot

fig = plt.figure(figsize=[16,12])

ax = fig.add_subplot(2,3,1,projection=amap)
amap.plot(clip_interval=[10,90]*u.percent)

ax = fig.add_subplot(2,3,4,projection=bmap)
bmap.plot(clip_interval=[10,90]*u.percent)

#ax = fig.add_subplot(2,3,3,projection=dmap)
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
pmap.plot(clip_interval=[20,90]*u.percent)

#======================================================================
# ratio

rat = np.where(bkg <= 0.0, 0.0, data/bkg)
rmap = sunpy.map.Map(rat,header)

emap = Enhance.fnrgf(rmap, instrument, detector)

ax = fig.add_subplot(2,3,6,projection=pmap)
pmap.plot(clip_interval=[20,90]*u.percent)

plt.show()