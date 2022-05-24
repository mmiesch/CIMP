"""
Checking the downsampling algorithm (reduce size by a factor of 2 in each dimension)
"""

import numpy as np
import matplotlib.pyplot as plt
import sunpy.map
import sunpy.visualization.colormaps as cm

from CIMP import Snapshot as snap

pcase = 1

# dcase = 1: subtract the background
# dcase = 2: take a ratio with the background
dcase = 1

clip = None
scales = None
rmask = None
colormap = 'stereo'

if pcase == 1:
    testcase = 1
    dcase = 2
    scale = (1,2)
elif pcase == 2:
    testcase = 2
    dcase = 2
    scale = (1,4)
else:
    print("specify a valid test case")
    exit()

#------------------------------------------------------------------------------

x = snap.snapshot.testcase(testcase)

# set dcase to 0 if there is no background
if dcase == 1:
    x.subtract_background()
elif dcase == 2:
    x.background_ratio(rescale = False)

map1 = x.map()

x.downsample(rescale = False)

map2 = x.map()

print(f"Size before downsample {map1.data.shape}")
print(f"Size after downsample {map2.data.shape}")


#------------------------------------------------------------------------------
# plot

if colormap == 'stereo':
    cmap = plt.get_cmap('stereocor2')
else:
    cmap = plt.get_cmap('soholasco2')

fig = plt.figure(figsize=[16,8])

print(f"image 1 range: {map1.min()} {map1.max()}")
ax = fig.add_subplot(1,2,1,projection=map1)
map1.plot(vmin=scale[0], vmax=scale[1],cmap=cmap)

print(f"image 2 range: {map2.min()} {map2.max()}")
ax = fig.add_subplot(1,2,2,projection=map2)
map2.plot(vmin=scale[0], vmax=scale[1],cmap=cmap)

plt.show()
