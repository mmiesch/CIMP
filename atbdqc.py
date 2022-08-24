
"""
Code for generating th QC figure in the CCOR ATBD
"""

import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import sunpy.visualization.colormaps as cm

from CIMP import Snapshot as snap
from CIMP import Enhance
from sunpy.net import attrs as a
#------------------------------------------------------------------------------

# L1 STEREO data
instrument1 = a.Instrument.secchi
detector1 = a.Detector.cor2
dir1='/home/mark.miesch/sunpy/data/secchi_cor2/L1/2012/09/'
bgfile1 = dir1+'background.fts'
file1 = dir1+'20/20120920_153900_14c2A.fts'
rmin = 0.16
rmax = 1.0

# L0.5 LASCO data
instrument2 = a.Instrument.lasco
detector2 = a.Detector.c3
dir2='/home/mark.miesch/data/lasco_monthly/c3/2012_04/'
bgfile2 = dir2+'background.fts'
file2 = dir2+'15/32296647.fts'
rmin = 0.16
rmax = 1.0

outdir = '/home/mark.miesch/Products/image_processing/figs/'

#------------------------------------------------------------------------------

fig = 1

if fig == 1:

    outfile = 'qc.png'

    cmap1 = plt.get_cmap('stereocor2')
    cmap2 = plt.get_cmap('soholasco2')

    title1 = 'Invalid pixels outside FOV'
    scale1 = (0, 1.0e-9)

    title2 = 'Corrupted image'
    scale2 = (0, 4000)

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

x1 = snap.snapshot(file = file1, bgfile = bgfile1, \
                  instrument = instrument1, detector = detector1)

#------------------------------------------------------------------------------
# second image

x2 = snap.snapshot(file = file2, bgfile = bgfile2, \
                  instrument = instrument2, detector = detector2)

#------------------------------------------------------------------------------
print(f"x1 minmax: {x1.min()} {x1.max()}")
print(f"x2 minmax: {x2.min()} {x2.max()}")

print(f"x1 res: {x1.data.shape[0]} {x1.data.shape[1]}")
print(f"x2 res: {x2.data.shape[0]} {x2.data.shape[1]}")

print(f"x1 time {x1.time}")
print(f"x2 time {x2.time}")

#------------------------------------------------------------------------------
# plot

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
