
"""
Code for generating some of the figures for the ATBD

"""

import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import sunpy.visualization.colormaps as cm

from CIMP import Snapshot as snap
from CIMP import Enhance
from sunpy.net import attrs as a
#------------------------------------------------------------------------------
"""
This replicates the process() method in the Animate class but with a little more freedom to specify the details of the processing

"""

def process(snap, background = 'ratio', point = 'none', \
            detail = 'none', noise = 'none', equalize = False, \
            downsample = False, clip = None, rmin = 0.0, rmax = None):

    if background == 'subtract':
        snap.mask_background(rmin = rmin, rmax = rmax, nonzero = False)
        snap.subtract_background(rescale=False)
    else:
        snap.mask_background(rmin = rmin, rmax = rmax, nonzero = True)
        snap.background_ratio(rescale=False)

    snap.mask_annulus(rmin = rmin, rmax = rmax)

    if downsample:
        snap.downsample()

    snap.enhance(clip = clip, point = point, detail = detail, noise_filter = noise, \
                 equalize = equalize)

    return

#------------------------------------------------------------------------------
# parameters common to all images

# L1 STEREO data
instrument = a.Instrument.secchi
detector = a.Detector.cor2
dir='/home/mark.miesch/sunpy/data/secchi_cor2/L1/2012/09/'
bgfile = dir+'background.fts'
file = dir+'20/20120920_153900_14c2A.fts'
rmin = 0.15
rmax = 1.0

outdir = '/home/mark.miesch/Products/image_processing/figs/'

#------------------------------------------------------------------------------
# choose the images you want to compare

fig = 3

if fig == 1:

    outfile = 'background_subtraction.png'

    cmap = plt.get_cmap('stereocor2')

    title1 = 'Background difference'
    background1 = 'subtract'
    downsample1 = False
    clip1       = None
    point1      = 'None'
    detail1     = 'None'
    noise1      = 'None'
    equalize1   = False
    scale1 = (0.0, 1.5e-10)

    title2 = 'Background ratio'
    background2 = 'ratio'
    downsample2 = False
    clip2       = None
    point2      = 'None'
    detail2     = 'None'
    noise2      = 'None'
    equalize2   = False
    scale2 = (1.0, 1.3)

elif fig == 2:

    outfile = 'downsampling.png'

    cmap = plt.get_cmap('stereocor2')

    title1 = 'Base image 2048 x 2048'
    background1 = 'ratio'
    downsample1 = False
    clip1       = None
    point1      = 'None'
    detail1     = 'None'
    noise1      = 'None'
    equalize1   = False
    scale1 = (1.0, 1.3)

    title2 = 'Downsample to 1024 x 1024'
    background2 = 'ratio'
    downsample2 = True
    clip2       = None
    point2      = 'None'
    detail2     = 'None'
    noise2      = 'None'
    equalize2   = False
    scale2 = (1.0, 1.3)

elif fig == 3:

    outfile = 'pf_median.png'

    cmap = plt.get_cmap('stereocor2')

    title1 = 'Base image'
    background1 = 'ratio'
    downsample1 = True
    clip1       = None
    point1      = 'None'
    detail1     = 'None'
    noise1      = 'None'
    equalize1   = False
    scale1 = (1.0, 1.3)

    title2 = 'Median bright point filter'
    background2 = 'ratio'
    downsample2 = True
    clip2       = None
    point2      = 'median'
    detail2     = 'None'
    noise2      = 'None'
    equalize2   = False
    scale2 = (1.0, 1.3)

else:
    print("pick a valid figure number")
    exit()

#------------------------------------------------------------------------------
# select color map

#cmap = plt.get_cmap('soholasco2')
#cmap = plt.get_cmap('stereocor2')

#------------------------------------------------------------------------------
# First image

x1 = snap.snapshot(file = file, bgfile = bgfile, \
                  instrument = instrument, detector = detector)

process(x1, background = background1, point = point1, detail = detail1, \
        noise = noise1, equalize = equalize1, downsample = downsample1, \
        clip = clip1, rmin = rmin, rmax = rmax)

if background1 == 'subtract':
    x1.powerlaw()

#------------------------------------------------------------------------------
# second image

x2 = snap.snapshot(file = file, bgfile = bgfile, \
                  instrument = instrument, detector = detector)

process(x2, background = background2, point = point2, detail = detail2, \
        noise = noise2, equalize = equalize2, downsample = downsample2, \
        clip = clip2, rmin = rmin, rmax = rmax)

#------------------------------------------------------------------------------
print(f"x1 minmax: {x1.min()} {x1.max()}")
print(f"x2 minmax: {x2.min()} {x2.max()}")

print(f"x1 res: {x1.data.shape[0]} {x1.data.shape[1]}")
print(f"x2 res: {x2.data.shape[0]} {x2.data.shape[1]}")

#------------------------------------------------------------------------------
# plot

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,6))

ax[0].imshow(x1.data,cmap=cmap,vmin=scale1[0],vmax=scale1[1], \
             origin='lower')
ax[0].axis('off')
ax[0].set_title(title1)

ax[1].imshow(x2.data,cmap=cmap,vmin=scale2[0],vmax=scale2[1], \
             origin='lower')
ax[1].axis('off')
ax[1].set_title(title2)

fig.tight_layout(pad=1,rect=(0.01,0.01,.99,.98))

plt.savefig(outdir+outfile)

plt.show()
