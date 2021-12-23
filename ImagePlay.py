"""
The purpose of this script is to start playing around with different python image processing tools
"""

import numpy as np
from CIMP import Event as ev
import sunpy.map
from sunpy.net import attrs as a
import matplotlib.pyplot as plt

from skimage import exposure

plotcase = 1

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

else:
    print("specify a valid plotcase")
    exit()    


x = ev.event.testcase(testcase)

# for experimenting
#timerange = a.Time('2016/09/06 8:00:00', '2016/09/06 12:00:00')
#x = ev.event.fromtime(a.Instrument.lasco, a.Detector.c2, timerange)

print(80*'-')
print(x)
print(80*'-')
print(repr(x))
print(80*'-')

# ===================
# pick an image to work with

a = x[idx]



# image as a sunpy map
amap = x.map(idx)

amin = np.amin(a)
amax = np.amax(a)

print(f"Original image range: {amin} to {amax}")

# clipped image
aclip = a.clip(min=scale[0], max=scale[1])

# scaled byte image
#asc = (a - amin)/(amax - amin)
asc = (255*aclip/np.amax(aclip)).astype('uint8')

# ---------
# basic enhancements

#psc = exposure.equalize_hist(asc)
psc = exposure.equalize_adapthist(asc)

p = (scale[1] - scale[0])*psc + scale[0]

#======================================================================
# plot

fig = plt.figure(figsize=[22,10])

if len(scale) != 2:
    scale = (amin, amax)

print(f"image scale: {scale[0]} to {scale[1]}")

# Optionally plot as a sunpy map

plotasmap = True

if plotasmap:
    ax = fig.add_subplot(1,2,1,projection=amap)
    amap.plot(vmin = scale[0], vmax = scale[1])
    
    pmap = sunpy.map.Map(p, x.header[0])
    ax = fig.add_subplot(1,2,2,projection=pmap)
    pmap.plot(vmin = scale[0], vmax = scale[1])

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