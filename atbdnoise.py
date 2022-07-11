
"""
This is similar to atbdfigs.py but focusing on noise removal.  Use a different set of images to emphasize the noise.

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
            downsample = False, clip = None, rmin = 0.0, rmax = None, \
            params = None):

    if background == 'subtract':
        snap.mask_background(rmin = rmin, rmax = rmax, nonzero = False)
        snap.subtract_background(rescale=False)
    elif background == 'ratio':
        snap.mask_background(rmin = rmin, rmax = rmax, nonzero = True)
        snap.background_ratio(rescale=False)
    else:
        snap.mask_background(rmin = rmin, rmax = rmax, nonzero = False)

    snap.mask_annulus(rmin = rmin, rmax = rmax)

    if downsample:
        snap.downsample()

    snap.enhance(clip = clip, point = point, detail = detail, noise_filter = noise, \
                 equalize = equalize, params = params)

    # hit it with another mask after processing
    snap.mask_annulus(rmin=rmin, rmax = rmax, missingval = np.nanmin(snap.data))

    return

#------------------------------------------------------------------------------
# parameters common to all images

# L0.5 LASCO data 
instrument = a.Instrument.lasco
detector = a.Detector.c3
dir='/home/mark.miesch/data/lasco_monthly/c3/2012_04/'
bgfile = dir+'background.fts'
file1 = dir+'15/32296620.fts'
file2 = file1
rmin = 0.16
rmax = 1.0
params1 = None
params2 = None

outdir = '/home/mark.miesch/Products/image_processing/figs/'

#------------------------------------------------------------------------------
# choose the images you want to compare

fig = 2

if fig == 1:

    outfile = 'pf_lasco_2012.png'

    cmap = plt.get_cmap('stereocor2')

    title1 = 'Median bright point filter'
    background1 = 'ratio'
    downsample1 = True
    clip1       = None
    point1      = 'median'
    detail1     = 'None'
    noise1      = 'None'
    equalize1   = False
    scale1 = (1.0, 1.3)

    title2 = 'OMR'
    background2 = 'ratio'
    downsample2 = True
    clip2       = None
    point2      = 'omr'
    detail2     = 'None'
    noise2      = 'None'
    equalize2   = False
    scale2 = (1.0, 1.3)

elif fig == 2:

    outfile = 'lasco_norm.png'

    cmap = plt.get_cmap('stereocor2')

    file1 = dir+'15/32296618.fts'
    title1 = file1.split('/')[-1]
    background1 = 'ratio'
    downsample1 = True
    clip1       = None
    point1      = 'None'
    detail1     = 'None'
    noise1      = 'None'
    equalize1   = False
    scale1 = (0.056, 0.074)

    file2 = dir+'15/32296619.fts'
    title2 = file2.split('/')[-1]
    background2 = 'ratio'
    downsample2 = True
    clip2       = None
    point2      = 'None'
    detail2     = 'None'
    noise2      = 'None'
    equalize2   = False
    scale2 = (0.056, 0.074)

else:
    print("pick a valid figure number")
    exit()

#------------------------------------------------------------------------------
# select color map

try:
    cmap1 = cmap
    cmap2 = cmap
except:
    pass

#------------------------------------------------------------------------------
# First image

x1 = snap.snapshot(file = file1, bgfile = bgfile, \
                  instrument = instrument, detector = detector)

process(x1, background = background1, point = point1, detail = detail1, \
        noise = noise1, equalize = equalize1, downsample = downsample1, \
        clip = clip1, rmin = rmin, rmax = rmax, params = params1)

print(f"file1 exposure time {x1.header['EXPTIME']}")

#------------------------------------------------------------------------------
# second image

x2 = snap.snapshot(file = file2, bgfile = bgfile, \
                  instrument = instrument, detector = detector)

process(x2, background = background2, point = point2, detail = detail2, \
        noise = noise2, equalize = equalize2, downsample = downsample2, \
        clip = clip2, rmin = rmin, rmax = rmax, params = params2)

print(f"file2 exposure time {x2.header['EXPTIME']}")

#------------------------------------------------------------------------------
print(f"x1 minmax: {x1.min()} {x1.max()}")
print(f"x2 minmax: {x2.min()} {x2.max()}")

print(f"x1 res: {x1.data.shape[0]} {x1.data.shape[1]}")
print(f"x2 res: {x2.data.shape[0]} {x2.data.shape[1]}")

print(f"x1 time {x1.header['DATE-OBS']}")

#------------------------------------------------------------------------------
# plot

if fig == 12:
    data1 = x1.data
    data2 = x1.background
else:
    data1 = x1.data
    data2 = x2.data

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,6))

ax[0].imshow(data1,cmap=cmap1,vmin=scale1[0],vmax=scale1[1], \
             origin='lower')
ax[0].axis('off')
ax[0].set_title(title1)

ax[1].imshow(data2,cmap=cmap2,vmin=scale2[0],vmax=scale2[1], \
             origin='lower')
ax[1].axis('off')
ax[1].set_title(title2)

fig.tight_layout(pad=1,rect=(0.01,0.01,.99,.98))

# label
label = True

if label:
    plt.annotate("(a)", (0.05,0.87), xycoords = 'figure fraction', color='white', \
                 fontsize = 'x-large', fontweight = 'semibold')

    plt.annotate("(b)", (0.54,0.87), xycoords = 'figure fraction', color='white', \
                 fontsize = 'x-large', fontweight = 'semibold')


plt.savefig(outdir+outfile)

plt.show()
