"""
A tool to plot out an event
"""

import numpy as np
from CIMP import Event as ev
from sunpy.net import attrs as a

#---------------------------------------------------

plotcase = 13

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
    scale = None

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
    scale = None
    clip = (0.0,200)

elif plotcase == 7:
    # Event 22 in NASA/NOAA MOU Annex Final Report (Mays et al 2015)
    testcase = 3
    instrument = a.Instrument.lasco
    detector = a.Detector.c2
    timerange = a.Time('2013/05/17 9:00:00', '2013/05/17 11:30:00')
    nrgf = False
    enhance = True
    plotframes = (1, 2, 4, 6)
    scale = (0, 100.0)
    clip = scale

elif plotcase == 8:
    # Event 22 in NASA/NOAA MOU Annex Final Report (Mays et al 2015)
    testcase = 3
    instrument = a.Instrument.lasco
    detector = a.Detector.c2
    timerange = a.Time('2013/05/17 9:00:00', '2013/05/17 11:30:00')
    nrgf = True
    enhance = True
    plotframes = (1, 2, 4, 6)
    scale = (0, 4.0)
    clip = scale

elif plotcase == 9:
    # Event 22 in NASA/NOAA MOU Annex Final Report (Mays et al 2015)
    testcase = 4
    instrument = a.Instrument.lasco
    detector = a.Detector.c3
    timerange = a.Time('2013/05/17 9:30:00', '2013/05/17 13:30:00')
    nrgf = False
    enhance = True
    plotframes = (4, 8, 12, 18)
    scale = None
    clip = (0,100.0)

elif plotcase == 10:
    # Event 22 in NASA/NOAA MOU Annex Final Report (Mays et al 2015)
    testcase = 5
    instrument = a.Instrument.secchi
    detector = a.Detector.cor1
    timerange = a.Time('2013/05/17 9:00:00', '2013/05/17 11:00:00')
    nrgf = False
    enhance = True
    plotframes = None
    scale = (0, 100.0)
    clip = scale

elif plotcase == 11:
    # Event 22 in NASA/NOAA MOU Annex Final Report (Mays et al 2015)
    testcase = 6
    instrument = a.Instrument.secchi
    detector = a.Detector.cor2
    timerange = a.Time('2013/05/17 9:30:00', '2013/05/17 13:00:00')
    nrgf = False
    enhance = False
    plotframes = (1,2,3,5)
    scale = (0, 300.0)
    clip = scale

elif plotcase == 12:
    # Event 22 in NASA/NOAA MOU Annex Final Report (Mays et al 2015)
    testcase = 7
    instrument = a.Instrument.secchi
    detector = a.Detector.cor1
    timerange = a.Time('2013/05/17 9:00:00', '2013/05/17 11:00:00')
    nrgf = False
    enhance = True
    plotframes = None
    scale = (0, 300.0)
    clip = scale

elif plotcase == 13:
    # Event 22 in NASA/NOAA MOU Annex Final Report (Mays et al 2015)
    testcase = 8
    instrument = a.Instrument.secchi
    detector = a.Detector.cor2
    timerange = a.Time('2013/05/17 9:30:00', '2013/05/17 13:00:00')
    nrgf = False
    enhance = True
    plotframes = (1,2,3,5)
    scale = None
    clip = (0,600)
else:
    print("specify a valid plotcase")
    exit()    

#---------------------------------------------------
# Create event

if testcase is None:

    # if testcase is not specified, you have to specify an instrument, detector, and time range
    x = ev.event.fromtime(instrument, detector, timerange)

else:
   x = ev.event.testcase(testcase)

print(80*'-')
print(x)
print(80*'-')
print(repr(x))
print(80*'-')

#---------------------------------------------------
# pick 4 frames to plot
if plotframes == None:
    n = np.uint16((x.nframes - 1)/4)
    plotframes = (n, 2*n, 3*n, x.nframes-1)

print(f"PLOTFRAMES: {plotframes}")

#---------------------------------------------------
# Optionally Apply a filter

if nrgf:
   x.nrgf()

if enhance:
    x.enhance(clip = clip, detail='mgn')
    #x.enhance(clip = clip, detail='contrast')
    #x.enhance(clip = clip, detail='fnrgf')
    #x.noise_gate(cubesize=(3,12,12)); scale=(0,100)

# ===================
import matplotlib.pyplot as plt
import astropy.units as u

fig = plt.figure(figsize=[18,10])

# ===================\

# plot the first frame, not as a diff image
amap = x.map(0)
ax = fig.add_subplot(2,3,1,projection=amap)
#amap.plot(clip_interval=[10,90]*u.percent)
amap.plot()

# pick a scale from one of the middle images
ref = x.map(plotframes[3])
print(f"Reference data range: {ref.min()} to {ref.max()}")

if scale is None:
    scale = (ref.min(), ref.max())

print(f"image scale: {scale[0]} to {scale[1]}")

for i in np.arange(0,4):
    amap = x.map(plotframes[i])
    ax = fig.add_subplot(2,3,i+2,projection=amap)
    if scale is None:
       pplot = amap.plot()
    else:
       pplot = amap.plot(vmin = scale[0], vmax = scale[1])
    print(f"color table {pplot.get_clim()}")

# plot sum in last frame
long_exposure = x.sum()

ax = fig.add_subplot(2,3,6,projection=long_exposure)
long_exposure.plot(clip_interval=[10,90]*u.percent)

plt.show()