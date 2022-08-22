
"""
This is similar to atbdfigs.py but focusing on noise removal.  Use a different set of images to emphasize the noise.

"""

import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import sunpy.visualization.colormaps as cm

from astropy.io import fits
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
            params = None, rescale_output = False):

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

    if rescale_output:
        snap.rescale()

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
rescale1 = False
rescale2 = False

outdir = '/home/mark.miesch/Products/image_processing/figs/'

#------------------------------------------------------------------------------
# choose the images you want to compare

fig = 10

x2file = None

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

elif fig == 3:

    outfile = 'lasco_noise_mgn.png'

    cmap = plt.get_cmap('soholasco2')

    file1 = dir+'15/32296635.fts'
    title1 = 'base image'
    background1 = 'ratio'
    downsample1 = False
    clip1       = None
    point1      = 'None'
    detail1     = 'None'
    noise1      = 'None'
    equalize1   = False
    scale1 = (1.0, 1.2)

    file2 = dir+'15/32296635.fts'
    title2 = "OMR/MGN enhanced"
    background2 = 'ratio'
    downsample2 = True
    clip2       = (1.0,1.4)
    point2      = 'omr'
    detail2     = 'mgn'
    params2     = (0.8,1.5)
    noise2      = 'None'
    equalize2   = False
    scale2 = (0.2,1.0)

elif fig == 4:

    outfile = 'lasco_noise_fnrgf.png'

    cmap = plt.get_cmap('stereocor2')

    file1 = dir+'15/32296635.fts'
    title1 = 'base image'
    background1 = 'ratio'
    downsample1 = False
    clip1       = None
    point1      = 'None'
    detail1     = 'None'
    noise1      = 'None'
    equalize1   = False
    scale1 = (1.0, 1.3)

    file2 = dir+'15/32296635.fts'
    title2 = "FNRGF enhanced"
    background2 = 'ratio'
    downsample2 = True
    clip2       = (1.0,1.3)
    point2      = 'None'
    detail2     = 'fnrgf'
    noise2      = 'None'
    equalize2   = False
    scale2 = (0.3666,0.371)

elif fig == 5:

    outfile = 'lasco_noise_mgn_omr.png'

    cmap = plt.get_cmap('stereocor2')

    file1 = dir+'15/32296635.fts'
    title1 = 'MGN enhanced'
    background1 = 'ratio'
    downsample1 = True
    clip1       = (1.0,1.3)
    point1      = 'None'
    detail1     = 'mgn'
    params1     = (0.8,1.2)
    noise1      = 'None'
    equalize1   = False
    scale1 = (0.2, 1.0)

    file2 = dir+'15/32296635.fts'
    title2 = "MGN enhanced with OMR"
    background2 = 'ratio'
    downsample2 = True
    clip2       = (1.0,1.3)
    point2      = 'omr'
    detail2     = 'mgn'
    params2     = (0.7,1.0)
    noise2      = 'None'
    equalize2   = False
    scale2 = (0.3,1.0)

elif fig == 6:

    outfile = 'lasco_noise_fnrgf.png'

    cmap = plt.get_cmap('stereocor2')

    file1 = dir+'15/32296635.fts'
    title1 = 'FNRGF enhanced'
    background1 = 'ratio'
    downsample1 = True
    clip1       = (1.0,1.3)
    point1      = 'None'
    detail1     = 'fnrgf'
    noise1      = 'None'
    equalize1   = False
    rescale1    = True
    scale1 = (0.1,1.0)

    file2 = dir+'15/32296635.fts'
    title2 = "FNRGF enhanced with OMR"
    background2 = 'ratio'
    downsample2 = True
    clip2       = (1.0,1.3)
    point2      = 'None'
    detail2     = 'fnrgf'
    noise2      = 'omr'
    equalize2   = False
    rescale2    = True
    scale2 = (0.5,1.0)

elif fig == 7:

    outfile = 'lasco_noise_bregman.png'

    cmap = plt.get_cmap('soholasco2')
    #cmap = plt.get_cmap('stereocor2')

    file1 = dir+'15/32296635.fts'
    title1 = 'OMR/MGN enhanced'
    background1 = 'ratio'
    downsample1 = True
    clip1       = (1.0,1.4)
    point1      = 'omr'
    detail1     = 'mgn'
    params1     = (0.8,1.5)
    noise1      = 'none'
    equalize1   = False
    scale1 = (0.2,1.0)

    file2 = dir+'15/32296635.fts'
    title2 = "with Bregman noise filter"
    background2 = 'ratio'
    downsample2 = True
    clip2       = (1.0,1.4)
    point2      = 'omr'
    detail2     = 'mgn'
    params2     = (0.8,1.2)
    noise2      = 'bregman'
    equalize2   = False
    scale2 = (0.1,1.0)

elif fig == 8:

    outfile = 'lasco_noise_omr.png'

    cmap = plt.get_cmap('soholasco2')
    #cmap = plt.get_cmap('stereocor2')

    file1 = dir+'15/32296635.fts'
    title1 = 'OMR/MGN enhanced'
    background1 = 'ratio'
    downsample1 = True
    clip1       = (1.0,1.4)
    point1      = 'omr'
    detail1     = 'mgn'
    params1     = (0.8,1.5)
    noise1      = 'none'
    equalize1   = False
    scale1 = (0.2,1.0)

    file2 = dir+'15/32296635.fts'
    title2 = "with OMR noise filter"
    background2 = 'ratio'
    downsample2 = True
    clip2       = (1.0,1.4)
    point2      = 'omr'
    detail2     = 'mgn'
    params2     = (0.8,1.2)
    noise2      = 'omr'
    equalize2   = False
    scale2 = (0.2,1.0)

elif fig == 9:

    outfile = 'lasco_noise_tv.png'

    cmap = plt.get_cmap('soholasco2')

    file1 = dir+'15/32296635.fts'
    title1 = 'OMR/MGN enhanced'
    background1 = 'ratio'
    downsample1 = True
    clip1       = (1.0,1.4)
    point1      = 'omr'
    detail1     = 'mgn'
    params1     = (0.8,1.5)
    noise1      = 'none'
    equalize1   = False
    scale1 = (0.2,1.0)

    file2 = dir+'15/32296635.fts'
    title2 = "with TV noise filter"
    background2 = 'ratio'
    downsample2 = True
    clip2       = (1.0,1.4)
    point2      = 'omr'
    detail2     = 'mgn'
    params2     = (0.8,1.2)
    noise2      = 'tv'
    equalize2   = False
    scale2 = (0.1,1.0)

if fig == 10:

    # The goal with this figure is to compare the
    # final product, after OMR, MGN, and NG, to the
    # original unprocessed image

    outfile = 'lasco_ng_pipeline.png'

    title1 = 'Base Image'
    file1 = dir+'15/32296635.fts'
    cmap1 = plt.get_cmap('soholasco2')
    background1 = 'ratio'
    downsample1 = False
    clip1       = None
    point1      = 'None'
    detail1     = 'None'
    noise1      = 'None'
    equalize1   = False
    scale1 = (1.0, 1.2)

    title2 = 'Processed image (OMR/MGN/NG)'
    cmap2 = plt.get_cmap('soholasco2')
    x2file = '/home/mark.miesch/data/lasco_monthly/c3/2012_04_ng/ng_hybrid_fig7.fts'
    scale2 = (0.2, 1.0)

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

x1 = snap.snapshot(file = file1, bgfile = bgfile, instrument = instrument, \
                   detector = detector, normalize = True)

process(x1, background = background1, point = point1, detail = detail1, \
        noise = noise1, equalize = equalize1, downsample = downsample1, \
        clip = clip1, rmin = rmin, rmax = rmax, params = params1, \
        rescale_output = rescale1)

print(f"file1 exposure time {x1.header['EXPTIME']}")

data1 = x1.data

#------------------------------------------------------------------------------
# second image

if x2file is None:

    x2 = snap.snapshot(file = file2, bgfile = bgfile, instrument = instrument, \
                       detector = detector, normalize = True)

    process(x2, background = background2, point = point2, detail = detail2, \
            noise = noise2, equalize = equalize2, downsample = downsample2, \
            clip = clip2, rmin = rmin, rmax = rmax, params = params2, \
            rescale_output = rescale2)

    print(f"file2 exposure time {x2.header['EXPTIME']}")

    data2 = x2.data

else:

    hdu = fits.open(x2file)[0]
    data2 = hdu.data

#------------------------------------------------------------------------------
print(f"x1 minmax: {np.min(data1)} {np.max(data1)}")
print(f"x2 minmax: {np.min(data2)} {np.max(data2)}")

print(f"x1 res: {data1.shape[0]} {data1.shape[1]}")
print(f"x2 res: {data2.shape[0]} {data2.shape[1]}")

print(f"x1 time {x1.time}")

#------------------------------------------------------------------------------
# plot

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
