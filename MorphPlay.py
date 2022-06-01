"""
Playing with Mathematical Morphology
"""

import numpy as np
import matplotlib.pyplot as plt
import sunpy.map
import sunpy.visualization.colormaps as cm

from CIMP import Snapshot as snap
from skimage.morphology import opening, closing, disk, erosion, reconstruction

pcase = 7

# dcase = 1: subtract the background
# dcase = 2: take a ratio with the background
dcase = 2

clip = None
scales = None
rmask = None
colormap = 'stereo'

if pcase == 1:
    comp = (0,1)
    testcase = 1
    scales = [(0,.05),(0, 1)]
elif pcase == 2:
    comp = (0,2)
    testcase = 1
    scales = [(0,.05),(0, 0.05)]
elif pcase == 3:
    comp = (0,3)
    testcase = 1
    scales = [(0,.05),(0, 0.05)]
elif pcase == 4:
    comp = (0,4)
    testcase = 1
    scales = [(0,.05),(0, 0.05)]
elif pcase == 5:
    comp = (0,5)
    testcase = 1
    scales = [(0,.05),(0, 0.05)]
elif pcase == 6:
    comp = (0,5)
    testcase = 2
    scales = [(0,.03),(0, 0.03)]
elif pcase == 7:
    comp = (0,6)
    testcase = 1
    scales = [(0,.05),(0.1, 1.0)]
    colormap = 'lasco'
else:
    print("specify a valid test case")
    exit()

#------------------------------------------------------------------------------
def plot_comparison(original, filtered, filter_name, scales, colormap = 'stereo'):
    if colormap == 'stereo':
        cmap = plt.get_cmap('stereocor2')
    else:
        cmap = plt.get_cmap('soholasco2')

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(16, 8), sharex=True,
                                   sharey=True)
    ax1.imshow(original, cmap=cmap, vmin = scales[0][0], vmax = scales[0][1])
    ax1.set_title('original')
    ax1.axis('off')
    ax2.imshow(filtered, cmap=cmap, vmin = scales[1][0], vmax = scales[1][1])
    ax2.set_title(filter_name)
    ax2.axis('off')

#------------------------------------------------------------------------------

x = snap.snapshot.testcase(testcase)

# set dcase to 0 if there is no background
if dcase == 1:
    x.subtract_background()
elif dcase == 2:
    x.background_ratio()

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

if comp.count(1) > 0:
    x.point_filter()
    titles.append("Point filter")
    images.append(x.data)
    dscales.append((0.0,1.0))

if comp.count(2) > 0:
    footprint = disk(3)
    fim = opening(x.data, footprint)
    titles.append("Opening footprint disk(3)")
    images.append(fim)
    dscales.append((0.0,1.0))

if comp.count(3) > 0:
    footprint = disk(3)
    fim = opening(x.data, footprint)
    fim2 = closing(fim, footprint)
    titles.append("Opening and closing")
    images.append(fim2)
    dscales.append((0.0,1.0))

if comp.count(4) > 0:
    footprint = disk(3)
    fim = erosion(x.data, footprint)
    titles.append("Erosion")
    images.append(fim)
    dscales.append((0.0,1.0))

if comp.count(5) > 0:
    footprint = disk(3)
    fim = erosion(x.data, footprint)
    fim2 = reconstruction(fim, x.data, method = 'dilation', footprint = footprint)
    titles.append("OMR")
    images.append(fim2)
    dscales.append((0.0,1.0))

if comp.count(6) > 0:
    #a = en.omr(x.data, rescaleim = False)
    x.enhance(point = 'omr', detail = 'mgn', noise_filter = 'omr')
    x.mask_annulus(rmax = 1.05)
    titles.append("Enhance")
    images.append(x.data)
    dscales.append((0.0,1.0))


#------------------------------------------------------------------------------

for image in images:
    print(f"range: {np.min(image)} {np.max(image)}")

if scales is None:
   scales = dscales
else:
    if scales[0] is None:
        scales[0] = dscales[0]
    if scales[1] is None:
        scales[1] = dscales[1]

#------------------------------------------------------------------------------
# plot

plot_comparison(images[0],images[1],titles[1], scales, colormap = colormap)

# ===================
# save to a file
# ==================

dir = '/home/mark.miesch/Products/image_processing/images/Morph/'

fname = f"Morph_t{testcase}_{comp[0]}_vs_{comp[1]}"

if dcase == 2:
    fname += "_rat"

if tag is not None:
    fname += f"_{tag}"

file = dir + fname + ".png"

plt.savefig(file)

plt.show()