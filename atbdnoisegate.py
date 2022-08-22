
"""
This script is specifically for making a noisegate figure for the ccor atbd
"""

import os
import matplotlib.pyplot as plt
import noisegate as ng
import numpy as np
import sunpy.visualization.colormaps as cm

from astropy.io import fits
from CIMP import Snapshot as snap
from CIMP import Enhance
from sunpy.net import attrs as a

# for warning / error statements; print red, yellow text to terminal
red = '\033[91m'
yellow = '\033[93m'
cend = '\033[0m'

#------------------------------------------------------------------------------
# parameters common to all images

# L0.5 LASCO data processed and qc filtered
dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04_qc'

# This is the image before noisegate is applied
file1 = dir+'/image068.fts'

outdir = '/home/mark.miesch/Products/image_processing/figs/'

#------------------------------------------------------------------------------

fig = 1

if fig == 1:

    outfile = 'lasco_ng_shot_f2.0.png'

    cmap = plt.get_cmap('soholasco2')
    #cmap = plt.get_cmap('stereocor2')

    title1 = 'OMR/MGN Enhanced'
    scale1 = (0.2, 1.0)

    title2 = 'Noisegate filter'
    cubesize = (12,12,12)
    model = 'shot'
    factor = 2.0
    dcidx = 67
    scale2 = (0.2, 1.0)

else:
    print("pick a valid figure number")
    exit()

#------------------------------------------------------------------------------
# First image - before filter

x1 = fits.open(file1)[0]

#------------------------------------------------------------------------------
# load data

# first pass: get list of all valid files
files = []
nx = None
ny = None
for file in os.listdir(dir):
    fpath = dir+'/'+file
    try:
        assert(os.path.isfile(fpath))
        assert("median" not in file)
        hdu = fits.open(fpath)[0]
        if nx is None:
            nx = hdu.data.shape[0]
            ny = hdu.data.shape[1]
        else:
            assert(hdu.data.shape[0] == nx)
            assert(hdu.data.shape[1] == ny)
        files.append(file)
    except Exception as e:
        print(red+f"{e}\nSkipping {file}"+cend)
        pass

# second pass: get data
Nfiles = len(files)
print(yellow+f"Number of valid files: {Nfiles}"+cend)
dcube = np.zeros((Nfiles,nx,ny))
idx = 0
for file in files:
    fpath = dir+'/'+file
    hdu = fits.open(fpath)[0]
    dcube[idx,:,:] = hdu.data.astype('float')
    idx += 1

print(yellow+f"Data shape: {dcube.shape}"+cend)

#------------------------------------------------------------------------------
# apply noisegate filter

ng.noise_gate_batch(dcube, cubesize=cubesize, model=model,
                    factor = factor)

#------------------------------------------------------------------------------
print(f"x1 time {x1.header['DATE-OBS']}")

#------------------------------------------------------------------------------
# plot

data1 = x1.data
data2 = dcube[dcidx,:,:]

print(f"x1 minmax: {np.min(data1)} {np.max(data1)}")
print(f"x2 minmax: {np.min(data2)} {np.max(data2)}")

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
