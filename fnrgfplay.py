
"""
This is for playing around with some of the parameters of the FNRGF to see what they do.

"""

import astropy.units as u
import sunkit_image.radial as radial
from skimage import exposure
from sunpy.map import Map
from sunkit_image.utils import equally_spaced_bins

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
def rescale(im):
    return exposure.rescale_intensity(im, out_range=(0,1))

def mask_annulus(im, rmin = 0.0, rmax = None, missingval = 0.0):
    """
    This sets the pixels inside rmin and/or outside rmax to the missing value (default 0)
    """

    nx = im.shape[0]
    ny = im.shape[1]

    nn = 0.5 * float(np.min((nx, ny)))

    rr1 = (rmin*nn)**2

    if rmax is None:
        rr2 = np.inf
    else:
        rr2 = (rmax*nn)**2

    x0 = 0.5*float(im.shape[0])
    y0 = 0.5*float(im.shape[1])

    for i in np.arange(0,im.shape[0]):
        for j in np.arange(0,im.shape[1]):
            r2 = (float(i)-x0)**2 + (float(j)-y0)**2
            if (r2 < rr1) or (r2 > rr2):
                im[i,j] = missingval


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

# L1 STEREO data
instrument = a.Instrument.secchi
detector = a.Detector.cor2
dir='/home/mark.miesch/sunpy/data/secchi_cor2/L1/2012/09/'
bgfile = dir+'background.fts'
file = dir+'20/20120920_153900_14c2A.fts'
#file = dir+'20/20120920_222400_14c2A.fts'
rmin = 0.16
rmax = 1.0
params1 = None
params2 = None

outdir = '/home/mark.miesch/Products/image_processing/figs/'

#------------------------------------------------------------------------------
# choose the images you want to compare

fig = 1

background = 'ratio'
downsample = True
clip       = (1.0,1.3)
point      = 'omr'
detail     = 'none'
noise      = 'none'
equalize   = False

if fig == 1:

    outfile = 'fnrgf1.png'

    cmap = plt.get_cmap('soholasco2')

    title1 = 'FNRGF1'
    scale1 = (0.15,1.0)
    kmax1 = 20
    rat1  = [1,15]
    n1 = 130
    c1 = radial.set_attenuation_coefficients(kmax1)

    title2 = 'FNRGF2'
    scale2 = (0.15,1.0)
    kmax2 = 10
    rat2  = [1,15]
    n2 = 130
    c2 = radial.set_attenuation_coefficients(kmax2)

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
# first do the downsampling and OMR point filter
# First image

x1 = snap.snapshot(file = file, bgfile = bgfile, \
                  instrument = instrument, detector = detector)

process(x1, background = background, point = point, detail = detail, \
        noise = noise, equalize = equalize, downsample = downsample, \
        clip = clip, rmin = rmin, rmax = rmax, params = params1)

#------------------------------------------------------------------------------
# parameters that should be the same for both

edges = equally_spaced_bins(3.0, 15.0) * u.R_sun

edges2 = equally_spaced_bins(1.0, 15.0) * u.R_sun

#------------------------------------------------------------------------------

input1 = x1.data
input2 = x1.data.copy()

header = x1.header

map1 = Map(input1, header)
Rmap1 = radial.fnrgf(map1, edges, kmax1, c1, ratio_mix = rat1)
data1 = rescale(Rmap1.data)
mask_annulus(data1, rmin=0.16, rmax = 1.0, missingval = np.nanmin(data1))

print(f"Map Scale: {map1.scale}")

map2 = Map(input2, header)
Rmap2 = radial.fnrgf(map2, edges, kmax2, c2, ratio_mix = rat2)
data2 = rescale(Rmap2.data)
mask_annulus(data2, rmin=0.16, rmax = 1.0, missingval = np.nanmin(data2))

#------------------------------------------------------------------------------
print(f"x1 minmax: {np.min(data1)} {np.max(data1)}")
print(f"x2 minmax: {np.min(data2)} {np.max(data2)}")

print(f"x1 res: {data1.shape[0]} {data2.shape[1]}")
print(f"x2 res: {data1.shape[0]} {data2.shape[1]}")

print(f"x1 time {header['DATE-OBS']}")

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
