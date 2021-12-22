"""
A tool to plot out an event
"""

import numpy as np
from CIMP import Event as ev
from sunpy.net import attrs as a

plotcase = 2

if plotcase == 1:
    testcase = 1
    nrgf = False
    plotframes = (3, 6, 9, 12)
    #scale = (700.0, 4000.0)
    scale = (0.0, 1000.0)

elif plotcase == 2:
    testcase = 1
    nrgf = True
    plotframes = (3, 6, 9, 12)
    scale = (0.0, 4.0)

else:
    print("specify a valid plotcase")
    exit()    

x = ev.event.testcase(1)

# for experimenting
#timerange = a.Time('2016/09/06 8:00:00', '2016/09/06 12:00:00')
#x = ev.event.fromtime(a.Instrument.lasco, a.Detector.c2, timerange)

print(80*'-')
print(x)
print(80*'-')
print(repr(x))
print(80*'-')

# ===================
# Optionally Apply a filter

if nrgf:
   x.nrgf()

# ===================
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u

fig = plt.figure(figsize=[18,10])

# plot the first frame, not as a diff image
amap = x.map(0)
ax = fig.add_subplot(2,3,1,projection=amap)
amap.plot(clip_interval=[10,90]*u.percent)

# pick 4 frames to plot
if len(plotframes) == 0:
    n = (x.nframes - 1)/4
    plotframes = [n, 2*n, 3*n, x.nframes-1] 

# pick a scale from one of the middle images
ref = x.map(9)
print(f"Reference data range: {ref.min()} to {ref.max()}")

if len(scale) != 2:
    scale = (ref.min(), ref.max())

print(f"image scale: {scale[0]} to {scale[1]}")

for i in np.arange(0,4):
    amap = x.map(plotframes[i])
    ax = fig.add_subplot(2,3,i+2,projection=amap)
    amap.plot(vmin = scale[0], vmax = scale[1])

# plot sum in last frame
long_exposure = x.sum()

ax = fig.add_subplot(2,3,6,projection=long_exposure)
long_exposure.plot(clip_interval=[10,90]*u.percent)

plt.show()