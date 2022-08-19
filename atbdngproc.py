
"""
The purpose of this script is to process a bunch of images for use with NoiseGate to model and remove noise.  By construction, it uses the same `process` function as `atbdnoise.py` but applies it to many images.
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import sunpy.visualization.colormaps as cm

from astropy.io import fits
from CIMP import Snapshot as snap
from CIMP import Enhance
from sunpy.net import attrs as a

#------------------------------------------------------------------------------

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
rmin = 0.16
rmax = 1.0
params = None
rescale = False

# pick the directories to do.  Note that the background is computed for
# the month so it is probably most accurate in the middle of the month
dlist = ['14','15','16']

outdir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04_processed/'

#------------------------------------------------------------------------------

fig = 1

if fig == 1:

    background = 'ratio'
    downsample = True
    clip       = (1.0,1.4)
    point      = 'omr'
    detail     = 'mgn'
    params     = (0.8,1.5)
    noise      = 'None'
    equalize   = False

else:
    print("pick a valid figure number")
    exit()

#------------------------------------------------------------------------------
# loop over all directories for the month

idx = np.uint64(0)
for d in dlist:
    day = dir+d
    for file in os.listdir(day):
        fpath = day+'/'+file
        try:
            assert(os.path.isfile(fpath))
            assert("median" not in file)
            assert(fpath != bgfile)
            x = snap.snapshot(file = fpath, bgfile = bgfile, \
                instrument = instrument, detector = detector, \
                normalize = True)
            process(x, background = background, point = point, detail = detail, \
                    noise = noise, equalize = equalize, downsample = downsample, \
                    clip = clip, rmin = rmin, rmax = rmax, params = params, \
                    rescale_output = rescale)
            outfile = outdir+'image'+str(idx).zfill(3)+'.fts'
            print(outfile)
            header0 = x.header
            del header0['HISTORY']
            hdu_out = fits.PrimaryHDU(x.data,header0)
            hdu_out.writeto(outfile, overwrite = True)
            idx += np.uint64(1)
        except:
            print(f"Skipping file {file}")
            pass

        hdu_out.writeto(outfile, overwrite = True)




#------------------------------------------------------------------------------
