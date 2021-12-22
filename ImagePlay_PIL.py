"""
The purpose of this script is to start playing around with different python image processing tools
"""

import numpy as np
from CIMP import Event as ev
import sunpy.map
from sunpy.net import attrs as a
import matplotlib.pyplot as plt

from PIL import Image
from PIL import ImageFilter
from PIL import ImageEnhance

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

# clipped array
#aclip = a.clip(min=0.0, max=1000.0)
aclip = a.clip(min=0.0, max=2000.0)

# image as a sunpy map
amap = x.map(idx)

# colored image as a numpy array:
#c = amap.cmap(a)

amin = np.amin(a)
amax = np.amax(a)

print(f"Original image range: {amin} to {amax}")

# ---------
# basic enhancements

#aimage = Image.fromarray(c.astype('uint8'),'RGB')
#aimage = Image.fromarray(np.uint8(c)).convert('RGB')

aim = Image.fromarray(a)

#aim = Image.fromarray(aclip)
#aim.show()
#argb = aim.convert('RGB')

agrey = aim.convert('L')
agrey.show()

#enh = ImageEnhance.Contrast(agrey)
#pim = enh.enhance(1.3)
#enh.enhance(1.3).show("more contrast")

enh = ImageEnhance.Sharpness(agrey)
pim = enh.enhance(2.0)

# processed image as numpy array
#pim = arargbgb.convert(mode = 'L')
pim.show()

# so this reproduces the original image
#pim = aim
#p = np.asarray(pim)

# this does not reproduce the original image
#pim = agrey
#parr = np.asarray(pim, dtype='float')/255
#p = parr*(amax - amin) + amin

#pim = agrey.convert(mode=None)
#p = np.asarray(pim)
parr = np.asarray(pim, dtype='float')/255
p = parr*(amax - amin) + amin

#plt.imshow(p,cmap=amap.cmap)
#plt.show()

print(f"Processed image range: {np.amin(p)} to {np.amax(p)}")

print(type(p))
print(p.shape)

# ---------
# plot


fig = plt.figure(figsize=[22,10])

if len(scale) != 2:
    scale = (amin, amax)

print(f"image scale: {scale[0]} to {scale[1]}")

ax = fig.add_subplot(1,2,1,projection=amap)
plt.imshow(a, vmin = scale[0], vmax = scale[1], cmap = amap.cmap)
print(f"Original image map: {amap.min()} to {amap.max()}")
#amap.plot(vmin = scale[0], vmax = scale[1])

# plot processed image

pmap = sunpy.map.Map(p, x.header[0])
print(f"Processed image map: {pmap.min()} to {pmap.max()}")
ax = fig.add_subplot(1,2,2,projection=pmap)
#pmap.plot(vmin = scale[0], vmax = scale[1])
#plt.imshow(p, vmin = scale[0], vmax = scale[1], cmap = amap.cmap)
plt.imshow(p,cmap=amap.cmap)


#plt.show()