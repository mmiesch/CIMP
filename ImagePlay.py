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
aclip = a.clip(min=0.0, max=1000.0)

# scaled byte image
#asc = (a - amin)/(amax - amin)
asc = (255*aclip/np.amax(aclip)).astype('uint8')

# ---------
# basic enhancements

#psc = exposure.equalize_hist(asc)
psc = exposure.equalize_adapthist(asc)

# ---------
# plot

fig = plt.figure(figsize=[22,10])

if len(scale) != 2:
    scale = (amin, amax)

print(f"image scale: {scale[0]} to {scale[1]}")

#ax = fig.add_subplot(1,2,2,projection=pmap)
#amap.plot(vmin = scale[0], vmax = scale[1])

ax = fig.add_subplot(1,2,1)
#plt.imshow(a, vmin = scale[0], vmax = scale[1], cmap = amap.cmap)
plt.imshow(asc, cmap=amap.cmap)

# plot processed image

#pmap = sunpy.map.Map(p, x.header[0])
#ax = fig.add_subplot(1,2,2,projection=pmap)
#pmap.plot(pmap)
#pmap.plot(vmin = scale[0], vmax = scale[1])
##plt.imshow(p, vmin = scale[0], vmax = scale[1], cmap = amap.cmap)
#plt.imshow(p,cmap=amap.cmap)

ax = fig.add_subplot(1,2,2)
plt.imshow(psc, cmap=amap.cmap)


plt.show()