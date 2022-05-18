"""
This is similar to EnhancePlay.py but it is a simplification, focusing on the development of the CIMP Snapshot class and the enhancement methods it calls through the CIMP Enhance class.  So, use EnhancePlay.py to assess new possible enhancement algorithms.  But, use this to test the algorithms that are already implemented in the Enhance class.
"""

import numpy as np
import matplotlib.pyplot as plt
import sunpy.map
import sunpy.visualization.colormaps as cm

from CIMP import Snapshot as snap

pcase = 4

# dcase = 1: subtract the background
# dcase = 2: take a ratio with the background
dcase = 1

scales = None
colormap = 'stereo'

if pcase == 1:
    comp = (0,1)
    testcase = 1
    dcase = 1
    pointfilter = True
    scales = [(0,1),(0.05,4.0)]
elif pcase == 2:
    comp = (0,2)
    testcase = 1
    dcase = 1
    pointfilter = True
    scales = [(0,1),(0.0,20.0)]
elif pcase == 3:
    comp = (0,3)
    testcase = 1
    dcase = 1
    pointfilter = True
    #scales = [(0,1),(-25,400.0)]
    colormap = 'lasco'
    scales = [(0.01,.12),(0,100.0)]
elif pcase == 4:
    comp = (0,4)
    testcase = 1
    dcase = 1
    pointfilter = False
    colormap = 'lasco'
    scales = [(0.01,.12),(0.0,1.0)]
elif pcase == 5:
    comp = (0,5)
    testcase = 1
    dcase = 1
    pointfilter = False
    colormap = 'lasco'
    scales = [(0.01,.12),(0.0,1.0)]
else:
    print("specify a valid test case")
    exit()

#------------------------------------------------------------------------------

x = snap.snapshot.testcase(testcase)

if dcase == 2:
    x.background_normalize()
else:
    x.subtract_background()

#------------------------------------------------------------------------------
# choose your battle

tag = None

images = []
titles = []

# default scales can be overridden
dscales = []

if comp.count(0) > 0:
    titles.append("Base image")
    images.append(x.data)
    dscales.append((0.0,1.0))

if pointfilter:
    x.point_filter()

if comp.count(1) > 0:
    titles.append("Point filter")
    images.append(x.data)
    dscales.append((0.0,1.0))

if comp.count(2) > 0:
    titles.append("NRGF")
    x.nrgf()
    images.append(x.data)
    dscales.append((0.0,4.0))

if comp.count(3) > 0:
    titles.append("FNRGF")
    x.fnrgf()
    images.append(x.data)
    dscales.append((0.0,4.0))

if comp.count(4) > 0:
    titles.append("enhance(mgn)")
    x.enhance(clip=(0.01, 0.12))
    x.mask_annulus(rmax = 1.05)
    images.append(x.data)
    dscales.append((0.0,1.0))

if comp.count(5) > 0:
    titles.append("enhance(fnrgf)")
    x.enhance(detail = 'fnrgf')
    images.append(x.data)
    dscales.append((0.0,1.0))

#------------------------------------------------------------------------------

if scales is None:
   scales = dscales
else:
    if scales[0] is None:
        scales[0] = dscales[0]
    if scales[1] is None:
        scales[1] = dscales[1]

#------------------------------------------------------------------------------
# plot

if colormap == 'stereo':
    cmap = plt.get_cmap('stereocor2')
else:
    cmap = plt.get_cmap('soholasco2')

fig = plt.figure(figsize=[16,8])

map1 = sunpy.map.Map(images[0],x.header)
print(f"image 1 range: {map1.min()} {map1.max()}")
ax = fig.add_subplot(1,2,1,projection=map1)
if scales[0] is None:
    map1.plot(title=titles[0],cmap=cmap)
else:
   print(f"image 1 scale: {scales[0][0]} to {scales[0][1]}")
   map1.plot(vmin=scales[0][0], vmax=scales[0][1], title=titles[0],cmap=cmap)

map2 = sunpy.map.Map(images[1],x.header)
print(f"image 2 range: {map2.min()} {map2.max()}")
ax = fig.add_subplot(1,2,2,projection=map2)
if scales[1] is None:
    map2.plot(title=titles[1],cmap=cmap)
else:
   print(f"image 2 scale: {scales[1][0]} to {scales[1][1]}")
   map2.plot(vmin=scales[1][0], vmax=scales[1][1], title=titles[1],cmap=cmap)

# ===================
# save to a file
# ===================
dir = '/home/mark.miesch/Products/image_processing/images/Enhance_snap/'

fname = f"Snap_t{testcase}_{comp[0]}_vs_{comp[1]}"

if dcase == 2:
    fname += "_rat"

if pointfilter:
    fname += "_pf"

if tag is not None:
    fname += f"_{tag}"

file = dir + fname + ".png"

plt.savefig(file)

plt.show()