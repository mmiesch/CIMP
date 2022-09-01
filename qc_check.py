"""
Simplified version of qcfilter_check.py that just loops over a bunch of image files,
applies a range filter, and flags them
"""
red = '\033[91m'
yellow = '\033[93m'
cend = '\033[0m'

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sunpy.visualization.colormaps as cm

from astropy.io import fits

def mask_annulus(im, rmin = 0.16, rmax = 1.0, missingval = 1.0):

    # annular radii in pixels (squared)
    nn = 0.5 * float(np.min((im.shape[0], im.shape[1])))
    rr1 = (rmin*nn)**2
    rr2 = (rmax*nn)**2

    # center of image
    x0 = 0.5*float(im.shape[0])
    y0 = 0.5*float(im.shape[1])

    x = np.tile(np.arange(im.shape[0]),(im.shape[1],1)).T
    y = np.tile(np.arange(im.shape[1]),(im.shape[0],1))
    rr = (x-x0)**2 + (y-y0)**2
    mask = (rr <= rr1) | (rr >= rr2)
    x = np.ma.masked_inside(rr, rr1, rr2)
    outside = float(x.count())/float(im.size)
    #print(f"outside FOV {outside}")
    return np.where(mask, missingval, im)

#------------------------------------------------------------------------------
def qc_range(image, range=(0.99,1.4), levels = (0.4, 0.5)):
    """
    Flag images outside of a given range
    range = desired range
    levels = threshold pixel fraction for flag values of 1 and 2
    """
    x = np.ma.masked_inside(image, range[0], range[1])
    outside = float(x.count())/float(image.size)
    if outside > levels[1]:
        flag = 2
    elif outside > levels[0]:
        flag = 1
    else:
        flag = 0

    return (flag, outside)

#------------------------------------------------------------------------------

fig = 1

if fig == 1:
    dir='/home/mark.miesch/Products/image_processing/ATBD/data/lasco_c3/L2proxy_2012_04'
elif fig == 2:

    dir='/home/mark.miesch/Products/image_processing/ATBD/data/lasco_c3/L2proxy_2014_01'

else:
    print("pick a valid figure number")
    exit()

#------------------------------------------------------------------------------

for file in os.listdir(dir):
    fpath = dir+'/'+file
    try:
        assert(os.path.isfile(fpath))
        hdu = fits.open(fpath)[0]
        x = mask_annulus(hdu.data)
        flag, invalid = qc_range(x)
        if flag > 0:
            print(80*'-'+f"\n{file}\n")
            print(f"   {flag} {invalid}")
    except Exception as e:
        print(red+f"{e}\nSkipping {file}"+cend)
        pass


