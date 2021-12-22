"""
The purpose of this script is to start playing around with different python image processing tools
"""

import numpy as np
from CIMP import Event as ev
import sunpy.map
from sunpy.net import attrs as a

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

# image as a sunpy map
amap = x.map(idx)

# colored image as a numpy array:
c = amap.cmap(a)

amin = np.amin(a)
amax = np.amax(a)

print(f"Original image range: {amin} to {amax}")

# ---------
# basic enhancements

#aimage = Image.fromarray(a.astype('uint8'),'RGB')
#aimage = Image.fromarray(np.uint8(a)).convert('RGB')

aim = Image.fromarray(a)
#argb = aim.convert('RGB')
#argb.show()
agrey = aim.convert('L')

#print(f"IMAGE MODES {aim.mode} {argb.mode}")
print(f"IMAGE MODES {aim.mode} {agrey.mode}")

#enh = ImageEnhance.Contrast(aim)
##pim = enh.enhance(1.3)
#enh.enhance(1.3).show("more contrast")

#enh = ImageEnhance.Sharpness(aim)
#asharp = enh.enhance(2.0)
#print(type(asharp))

# processed image as numpy array
#pim = aim
#pim = arargbgb.convert(mode = 'L')
pim = agrey
pim.show()
parr = np.asarray(pim, dtype='float')/255
p = parr*(amax - amin) + amin

print(f"Processed image range: {np.amin(p)} to {np.amax(p)}")

print(type(p))
print(p.shape)

# ---------
# plot

import matplotlib.pyplot as plt

fig = plt.figure(figsize=[22,10])

if len(scale) != 2:
    scale = (amin, amax)

print(f"image scale: {scale[0]} to {scale[1]}")

ax = fig.add_subplot(1,2,1,projection=amap)
#plt.imshow(a, vmin = scale[0], vmax = scale[1])
print(f"Original image map: {amap.min()} to {amap.max()}")
amap.plot(vmin = scale[0], vmax = scale[1])

# plot processed image

pmap = sunpy.map.Map(p, x.header[0])
print(f"Processed image map: {pmap.min()} to {pmap.max()}")
ax = fig.add_subplot(1,2,2,projection=pmap)
pmap.plot(vmin = scale[0], vmax = scale[1])

plt.show()