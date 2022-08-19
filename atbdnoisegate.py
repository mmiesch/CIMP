
"""
This script is specifically for making a noisegate figure for the ccor atbd
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
# parameters common to all images

# L0.5 LASCO data processed and qc filtered
dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04_qc'

# This is the image before noisegate is applied
file1 = dir+'/image068.fts'

outdir = '/home/mark.miesch/Products/image_processing/figs/'

#------------------------------------------------------------------------------
# choose the images you want to compare

fig = 1

if fig == 1:

    outfile = 'lasco_ng1.png'

    cmap = plt.get_cmap('soholasco2')
    #cmap = plt.get_cmap('stereocor2')

    title1 = 'OMR/MGN Enhanced'
    scale1 = (0.2, 1.0)

    title2 = 'Noisegate filter'
    scale2 = (0.2, 1.0)

else:
    print("pick a valid figure number")
    exit()

#------------------------------------------------------------------------------
# First image

x1 = fits.open(file1)[0]

#------------------------------------------------------------------------------
# second image

#------------------------------------------------------------------------------
print(f"x1 minmax: {np.min(x1.data)} {np.max(x1.data)}")
#print(f"x2 minmax: {x2.min()} {x2.max()}")

print(f"x1 res: {x1.data.shape[0]} {x1.data.shape[1]}")
#print(f"x2 res: {x2.data.shape[0]} {x2.data.shape[1]}")

print(f"x1 time {x1.header['DATE-OBS']}")

#------------------------------------------------------------------------------
# plot

data1 = x1.data
data2 = x1.data
#data2 = x2.data

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,6))

ax[0].imshow(data1,cmap=cmap,vmin=scale1[0],vmax=scale1[1], \
             origin='lower')
ax[0].axis('off')
ax[0].set_title(title1)

ax[1].imshow(data2,cmap=cmap,vmin=scale2[0],vmax=scale2[1], \
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
