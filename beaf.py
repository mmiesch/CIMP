
"""
Script to create movies showing images before and after processing
"""

import numpy as np
import os

from astropy.io import fits

#------------------------------------------------------------------------------
rootdir = '/home/mark.miesch/Products/image_processing/ATBD'
noisegate = True
framedir = None

#------------------------------------------------------------------------------
# get a few images just for testing the plot layout

dir = rootdir + '/data/lasco_c3/L3_2021_05'

files = []
files.append('LASCOC3_L3_2021_05_17_001807.fts')
files.append('LASCOC3_L3_2021_05_17_003007.fts')

hdu = fits.open(dir+'/'+files[0])
im1 = hdu[0].data
hdu.close()

hdu = fits.open(dir+'/'+files[1])
im2 = hdu[1].data
hdu.close()

#------------------------------------------------------------------------------
