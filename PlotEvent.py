"""
A tool to plot out an event
"""

import numpy as np
from CIMP import Event as ev
from sunpy.net import attrs as a

x = ev.event.testcase(1)

#timerange = a.Time('2016/09/06 8:00:00', '2016/09/06 12:00:00')
#x = ev.event.fromtime(a.Instrument.lasco, a.Detector.c2, timerange)

print(80*'-')
print(x)
print(80*'-')
print(repr(x))
print(80*'-')


# ===================
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u

fig = plt.figure(figsize=[18,8])

# plot the first frame, not as a diff image
amap = x.map(0)
ax = fig.add_subplot(2,3,1,projection=amap)
amap.plot(clip_interval=[10,90]*u.percent)

# pick 4 frames to plot 
p = [3,6,9,12]

# pick a scale from one of the middle images
#ref = x.map(9)
#vmin = ref.min()
#vmax = ref.max()
vmin = 700.0
vmax = 4000.0

print(f"image scale: {vmin} to {vmax}")

for i in np.arange(0,4):
    amap = x.map(p[i])
    ax = fig.add_subplot(2,3,i+2,projection=amap)
    amap.plot(vmin = vmin, vmax = vmax)

# plot sum in last frame
long_exposure = x.sum()

ax = fig.add_subplot(2,3,6,projection=long_exposure)
long_exposure.plot(clip_interval=[10,90]*u.percent)

plt.show()