
"""
Script to create movies showing images before and after processing
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import sunpy.visualization.colormaps as cm

from astropy.io import fits

#------------------------------------------------------------------------------
fig = 1

rootdir = '/home/mark.miesch/Products/image_processing/ATBD'
ngflag = False
framedir = None

if fig == 1:
    source = 'lasco'
    dir2 = rootdir + '/data/lasco_c3/L3_2021_05'
    cmap2 = plt.get_cmap('soholasco2')
    endfile2 = 'LASCOC3_L3_2021_05_17_013020.fts'
    duration = 2.0
    ngflag = False
    scale2 = (0.0, 0.5)

    dir1 = rootdir + '/data/lasco_c3/L2proxy_2021_05'
    cmap1 = plt.get_cmap('stereocor2')
    scale1 = (1.0, 1.3)

    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/lasco_2021_05_16_ba.mp4'
    framedir = pdir+f'/frames/debug'

else:
    print("pick a valid figure number")
    exit()

#------------------------------------------------------------------------------

#fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,6))
#
#ax[0].imshow(im1,cmap=cmap1,vmin=scale1[0],vmax=scale1[1], \
#             origin='lower')
#ax[0].axis('off')
#
#ax[1].imshow(im2,cmap=cmap2,vmin=scale2[0],vmax=scale2[1], \
#             origin='lower')
#ax[1].axis('off')
#
#fig.tight_layout(pad=1,rect=(0.01,0.01,.99,.99))

#plt.show()