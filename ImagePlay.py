"""
The purpose of this script is to start playing around with different python image processing tools
"""

import numpy as np
from CIMP import Event as ev
from CIMP import Snapshot as snap
import sunpy.map
from sunpy.net import attrs as a
import matplotlib.pyplot as plt

from skimage import exposure
from skimage.filters.rank import (median, enhance_contrast, enhance_contrast_percentile,
                                  autolevel)
from skimage.morphology import disk
from skimage.restoration import (denoise_tv_chambolle, denoise_nl_means)

plotcase = 3

snapshot = False

if plotcase == 1:
    testcase = 1
    nrgf = False
    idx = 9
    scale = (0.0, 1000.0)

elif plotcase == 2:
    testcase = 1
    nrgf = True
    idx = 9
    scale = (0.0, 4.0)

elif plotcase == 3:
    snapshot = True
    testcase = 1
    nrgf = False
    scale = (600, 4000)

else:
    print("specify a valid plotcase")
    exit()    

if snapshot:
    x = snap.snapshot.testcase(testcase)
    a = x.data
    amap = x.map()
    header = x.header
else:
    x = ev.event.testcase(testcase)
    a = x[idx]
    amap = x.map(idx)
    header = x.header[idx]

    # for experimenting
    #timerange = a.Time('2016/09/06 8:00:00', '2016/09/06 12:00:00')
    #x = ev.event.fromtime(a.Instrument.lasco, a.Detector.c2, timerange)

print(80*'-')
print(x)
print(80*'-')
print(repr(x))
print(80*'-')

# ===================

amin = np.amin(a)
amax = np.amax(a)

if len(scale) != 2:
    scale = (amin, amax)

print(f"Original image range: {amin} to {amax}")

# clipped image
aclip = a.clip(min=scale[0], max=scale[1])

# scaled byte image
#asc = (a - amin)/(amax - amin)
#asc = (255*aclip/np.amax(aclip)).astype('uint8')
asc = aclip/np.amax(aclip)

#======================================================================
# image enhancements

# first remove noise
psc1 = denoise_tv_chambolle(asc, weight = 0.2)

# this doesn't do much
#psc = exposure.equalize_hist(asc)

# adaptive equalization is great
psc2 = exposure.equalize_adapthist(psc1)

# gamma correction
#psc = exposure.adjust_gamma(psc2, gamma=0.9, gain=1.0)

# enhance contrast
psc = enhance_contrast_percentile(psc2,disk(1), p0=.1, p1=.9)
#psc = autolevel(psc2,disk(10))

# rescale
psc = (psc.astype('float') - np.amin(psc))/(np.amax(psc) - np.amin(psc))

p = (scale[1] - scale[0])*psc + scale[0]

#======================================================================
# plot

fig = plt.figure(figsize=[22,10])

print(f"image scale: {scale[0]} to {scale[1]}")

# Optionally plot as a sunpy map

plotasmap = True

if plotasmap:
    ax = fig.add_subplot(1,2,1,projection=amap)
    amap.plot(vmin = scale[0], vmax = scale[1])
    
    pmap = sunpy.map.Map(p, header)
    ax = fig.add_subplot(1,2,2,projection=pmap)
    pplot = pmap.plot(vmin = scale[0], vmax = scale[1])

else:
    ax = fig.add_subplot(1,2,1)
    plt.imshow(a, vmin = scale[0], vmax = scale[1], cmap = amap.cmap)
    #plt.imshow(asc, cmap=amap.cmap)
    
    # plot processed image
    ax = fig.add_subplot(1,2,2)
    pplot = plt.imshow(psc, cmap=amap.cmap)
    
    # see what matplotlib is using for the limits of the color scale
    
clim = pplot.get_clim()
    
print(f"Limits used for the color table = {clim[0]} to {clim[1]}")

#======================================================================
plt.show()