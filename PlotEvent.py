"""
A tool to plot out an event
"""

import numpy as np
from CIMP import Event as ev
from sunpy.net import attrs as a

plotcase = 6

if plotcase == 1:
    testcase = 1
    nrgf = False
    enhance = False
    plotframes = (3, 6, 9, 12)
    scale = (0.0, 1000.0)

elif plotcase == 2:
    testcase = 1
    nrgf = True
    enhance = False
    plotframes = (3, 6, 9, 12)
    scale = (0.0, 4.0)

elif plotcase == 3:
    testcase = 1
    nrgf = False
    enhance = True
    plotframes = (3, 6, 9, 12)
    clip = (0.0, 1000.0)
    scale = (0.0, 1000.0)

elif plotcase == 4:
    testcase = 1
    nrgf = True
    enhance = True
    plotframes = (3, 6, 9, 12)
    scale = (0.0, 4.0)
    clip = scale

elif plotcase == 5:
    testcase = 2
    nrgf = True
    enhance = True
    plotframes = (7, 14, 21, 29)
    scale = (0.0, 4.0)
    clip = scale

elif plotcase == 6:
    testcase = 2
    nrgf = False
    enhance = True
    plotframes = (7, 14, 21, 29)
    scale = (0.0, 100.0)
    clip = scale

else:
    print("specify a valid plotcase")
    exit()    

x = ev.event.testcase(testcase)

# for experimenting
#timerange = a.Time('2016/09/06 8:00:00', '2016/09/06 12:00:00')
#timerange = a.Time('2016/09/06 10:00:00', '2016/09/06 18:30:00')
#x = ev.event.fromtime(a.Instrument.lasco, a.Detector.c3, timerange)
#exit()

print(80*'-')
print(x)
print(80*'-')
print(repr(x))
print(80*'-')

# pick 4 frames to plot
if plotframes == None:
    n = np.uint16((x.nframes - 1)/4)
    plotframes = (n, 2*n, 3*n, x.nframes-1)

print(f"PLOTFRAMES: {plotframes}")

# ===================
# Optionally Apply a filter

if nrgf:
   x.nrgf()

if enhance:
    x.enhance(clip = clip)

# ===================
import matplotlib.pyplot as plt
import astropy.units as u

fig = plt.figure(figsize=[18,10])

# ===================\

# plot the first frame, not as a diff image
amap = x.map(0)
ax = fig.add_subplot(2,3,1,projection=amap)
amap.plot(clip_interval=[10,90]*u.percent)

# pick a scale from one of the middle images
ref = x.map(plotframes[3])
print(f"Reference data range: {ref.min()} to {ref.max()}")

if len(scale) != 2:
    scale = (ref.min(), ref.max())

print(f"image scale: {scale[0]} to {scale[1]}")

for i in np.arange(0,4):
    print(f"MSM {type(i)} {type(plotframes[i])}")
    amap = x.map(plotframes[i])
    ax = fig.add_subplot(2,3,i+2,projection=amap)
    pplot = amap.plot(vmin = scale[0], vmax = scale[1])
    print(f"color table {pplot.get_clim()}")

# plot sum in last frame
long_exposure = x.sum()

ax = fig.add_subplot(2,3,6,projection=long_exposure)
long_exposure.plot(clip_interval=[10,90]*u.percent)

plt.show()