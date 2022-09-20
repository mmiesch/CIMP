
"""
The purpose of this script is to time how long the
noisegate processing takes.
As a baseline, consider a 1-week movie with the cadence
of CCOR (15 min).  This corresponds to 672 frames.
"""
import os
import numpy as np
import noisegate as ng

from astropy.io import fits
from time import perf_counter

#------------------------------------------------------------------------------

dir = '/home/mark.miesch/Products/image_processing/ATBD/data/timings'

#Nfiles = 96
Nfiles = 672

#------------------------------------------------------------------------------

dirlist = os.listdir(dir)
flist = list(sorted(dirlist, reverse=True))
idx = 0

files = []
while (len(files) < Nfiles):
    hdu = fits.open(dir+'/'+flist[idx])
    try:
        flag = hdu[0].header['L3QCFLAG']
    except:
        flag = 0
    if flag < 2:
        print(f"{flist[idx]}")
        files.append(flist[idx])
    else:
        print(f"file skipped {flist[idx]}")
    idx += 1

assert Nfiles == len(files)

#------------------------------------------------------------------------------

hdu = fits.open(dir+'/'+files[0])
nx, ny = hdu[0].data.shape
hdu.close()

images = np.zeros((Nfiles,nx,ny), dtype = 'float')
for idx in np.arange(Nfiles):
   hdu = fits.open(dir+'/'+files[idx])
   images[Nfiles-1-idx,:,:] = hdu[0].data
   hdu.close()


tstart = perf_counter()

out = ng.noise_gate_batch(images, cubesize=(18,18,18), model='constant', \
                       factor = 6.0)

tstop = perf_counter()

dt = (tstop - tstart) / 60.0

print(f"Total Time (min)   : {dt}")

#------------------------------------------------------------------------------
