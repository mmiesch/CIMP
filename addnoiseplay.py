"""
Add noise to CME model images
"""

import numpy as np
import matplotlib.pyplot as plt
import sunpy.visualization.colormaps as cm

from astropy.io import fits
from skimage.util import random_noise

dir = '/home/mark.miesch/Products/image_processing/ATBD/data/model/CME0_pos30'

indir = dir+'/L2proxy'

#cmap = plt.get_cmap('soholasco2')
cmap = plt.get_cmap('stereocor2')

#------------------------------------------------------------------------------
# read input image

infile = indir+'/Model0_2010_04_14_150634.fts'

clip = (1.0,4.0)

hdu = fits.open(infile)
data = hdu[0].data
hdu.close()
data = (data - clip[0])/(clip[1] - clip[0])
data = data.clip(min = 0.0, max = 1.0)

#------------------------------------------------------------------------------
ccase = 1

if ccase == 1:
  noisyd = random_noise(data, mode = 'gaussian', var = 0.05)
elif ccase == 2:
  noisyd = random_noise(data, mode = 'poisson')
else:
    print("pick a valid ccase")
    exit()

#------------------------------------------------------------------------------
# plot

print(f"{np.min(noisyd)} {np.max(noisyd)}")

plt.imshow(noisyd, vmin=0.0, vmax=1.0, cmap = cmap)

plt.show()
