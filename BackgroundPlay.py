"""
Experiment with different types of background images and see how they affect the generation of "pretties"
"""

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from sunpy.net import attrs as a
from sunpy.net import Fido
import sunpy.io
import sunpy.map

# basic image to work with
#file = '/home/mark.miesch/data/lasco_ssw/32352255.fts'
#file = '/home/mark.miesch/sunpy/data/lasco_c3/32305540.fts'
file = '/home/mark.miesch/sunpy/data/lasco_c3/32305543.fts'
data, header = sunpy.io.fits.read(file)[0]
amap = sunpy.map.Map(data,header)
print(f"data range: {amap.min()}, {amap.max()}")

# Background file provided by NRL
#bfile = '/home/mark.miesch/data/lasco_ssw/3m_clcl_130910.fts'
bfile = '/home/mark.miesch/data/lasco_ssw/3m_clcl_120716.fts'
bkg, bheader = sunpy.io.fits.read(bfile)[0]
#bmap = sunpy.map.Map(bkg,bheader)
bmap = sunpy.map.Map(bkg,header)
print(f"bkg range: {bmap.min()}, {bmap.max()}")

# use earlier image as a background
#dfile = '/home/mark.miesch/data/lasco_ssw/32352249.fts'
#dfile = '/home/mark.miesch/sunpy/data/lasco_c3/32305537.fts'
dfile = '/home/mark.miesch/sunpy/data/lasco_c3/32305537.fts'
dimg, dheader = sunpy.io.fits.read(dfile)[0]
dmap = sunpy.map.Map(dimg,dheader)
print(f"diff range: {dmap.min()}, {dmap.max()}")

# Now plot

fig = plt.figure(figsize=[22,10])

ax = fig.add_subplot(2,3,1,projection=amap)
amap.plot(clip_interval=[10,90]*u.percent)

ax = fig.add_subplot(2,3,4,projection=bmap)
bmap.plot(clip_interval=[10,90]*u.percent)

#ax = fig.add_subplot(2,3,3,projection=dmap)
#dmap.plot(clip_interval=[10,90]*u.percent)

#======================================================================
# differences in second column

pa = data - dimg
pmap = sunpy.map.Map(pa,header)
ax = fig.add_subplot(2,3,2,projection=pmap)
pmap.plot(clip_interval=[10,90]*u.percent)

pa = data - bkg
pmap = sunpy.map.Map(pa,header)
ax = fig.add_subplot(2,3,5,projection=pmap)
pmap.plot(clip_interval=[10,90]*u.percent)

#======================================================================
# ratio in third column

rat = np.where(dimg <= 0.0, 0.0, data/dimg)
rmap = sunpy.map.Map(rat,header)
ax = fig.add_subplot(2,3,3,projection=rmap)
rmap.plot(clip_interval=[10,90]*u.percent)

rat = np.where(bkg <= 0.0, 0.0, data/bkg)
rmap = sunpy.map.Map(rat,header)
ax = fig.add_subplot(2,3,6,projection=rmap)
rmap.plot(clip_interval=[0,70]*u.percent)

plt.show()