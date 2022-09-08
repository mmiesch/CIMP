
"""
Script to create movies showing images before and after processing
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import sunpy.visualization.colormaps as cm

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
im2 = hdu[0].data
hdu.close()

cmap1 = plt.get_cmap('soholasco2')
scale1 = (0.0,1.0)

cmap2 = plt.get_cmap('soholasco2')
scale2 = (0.0,1.0)

#------------------------------------------------------------------------------

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,6))

ax[0].imshow(im1,cmap=cmap1,vmin=scale1[0],vmax=scale1[1], \
             origin='lower')
ax[0].axis('off')

ax[1].imshow(im2,cmap=cmap2,vmin=scale2[0],vmax=scale2[1], \
             origin='lower')
ax[1].axis('off')

fig.tight_layout(pad=1,rect=(0.01,0.01,.99,.99))

plt.show()