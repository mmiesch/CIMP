
"""
Code for generating some of the figures for the ATBD
This just takes two fits files and plots the images
"""

import matplotlib.pyplot as plt
import numpy as np
import sunpy.visualization.colormaps as cm

from astropy.io import fits
from CIMP.Enhance import mask_annulus

#------------------------------------------------------------------------------
# choose the images you want to compare

outdir = '/home/mark.miesch/Products/image_processing/figs/'

fig = 1

if fig == 1:

    outfile = 'qclevels.png'

    title1 = 'QC flag = 1'
    dir1='/home/mark.miesch/Products/image_processing/ATBD/data/lasco_c3/L2proxyb_2014_01'
    file1 = dir1+'/'+ 'LASCOC3_2014_01_16_233006.fts'
    cmap1 = plt.get_cmap('stereocor2')
    #cmap1 = plt.get_cmap('soholasco2')
    scale1 = (0.98, 1.3)

    title2 = 'QC flag = 2'
    dir2 = dir1
    file2 = dir2+'/'+ 'LASCOC3_2014_01_13_185406.fts'
    cmap2 = plt.get_cmap('soholasco2')
    scale2 = (0.98, 1.3)

else:
    print("pick a valid figure number")
    exit()

#------------------------------------------------------------------------------
# load images

hdu1 = fits.open(file1)
data1 = hdu1[0].data
hdu1.close()

hdu2 = fits.open(file2)
data2 = hdu2[0].data
hdu2.close()

#------------------------------------------------------------------------------
# mask annulus

data1 = mask_annulus(data1,rmax = 1.0)
data2 = mask_annulus(data2,rmax = 1.0)

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
