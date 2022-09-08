
"""
Script to create movies showing images before and after processing
"""

import datetime
import matplotlib.pyplot as plt
import noisegate as ng
import numpy as np
import os
import sunpy.visualization.colormaps as cm

from astropy.io import fits

#------------------------------------------------------------------------------
def get_time(header, source):

    if source == 'lasco':
        d = np.array(header['DATE-OBS'].split('/'),dtype='int')
        t = np.rint(np.array(header['TIME-OBS'].split(':'),dtype='float')).astype('int')
        return datetime.datetime(d[0],d[1],d[2],t[0],t[1],t[2])
    elif source == 'stereo':
        d = header['DATE-OBS']
        return datetime.datetime.fromisoformat(d)
    else:
        print("ERROR: cannot find time")
        return 0

#------------------------------------------------------------------------------
def noisegate(images, cubesize=(12,12,12), factor = 4.0):

    nap = int((2*cubesize[0])/3)
    nt, nx, ny = images.shape
    dcube = np.zeros((nt+nap, nx, ny), dtype = 'float')
    dcube[0:nt,:,:] = images
    dcube[nt:,:,:] = np.flip(images[nt-nap:,:,:], axis = 0)

    print(f"Applying noisegate {images.shape} {dcube.shape}")

    dcubeng = ng.noise_gate_batch(dcube, cubesize=cubesize, model='constant',
                factor = factor)

    return dcubeng[nap:-nap,:,:]

#------------------------------------------------------------------------------
fig = 1

rootdir = '/home/mark.miesch/Products/image_processing/ATBD'
ngflag = False
framedir = None

if fig == 1:
    source = 'lasco'
    dir2 = rootdir + '/data/lasco_c3/L3_2021_05'
    cmap2 = plt.get_cmap('soholasco2')
    endfile = 'LASCOC3_L3_2021_05_17_013020.fts'
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
# Compile list of valid L3 files in time range of interest,
# with time stamps

dirlist = os.listdir(dir2)
flist = list(sorted(dirlist, reverse=True))
idx = flist.index(endfile)

dtmax = datetime.timedelta(days=duration)
dt = datetime.timedelta(days=0.0)

files = []
times = []
while (dt <= dtmax) and (idx < len(flist)):
    hdu = fits.open(dir2+'/'+flist[idx])
    try:
        flag = hdu[0].header['L3QCFLAG']
    except:
        flag = 0
    if flag < 2:
        t = get_time(hdu[0].header, source)
        if len(times) > 0:
            dt = times[0] - t
        print(f"{flist[idx]} {dt}")
        files.append(flist[idx])
        times.append(t)
    idx += 1
#------------------------------------------------------------------------------
# load images and apply noisegate

Nfiles = len(files)

hdu2 = fits.open(dir2+'/'+files[0])
nx2, ny2 = hdu2[0].data.shape
hdu2.close()

file1 = files[0].replace('_L3','')
hdu1 = fits.open(dir1+'/'+file1)
nx1, ny1 = hdu1[0].data.shape
hdu1.close()

images1 = np.zeros((Nfiles,nx1,ny1), dtype = 'float')
images2 = np.zeros((Nfiles,nx2,ny2), dtype = 'float')
for idx in np.arange(Nfiles):
   hdu2 = fits.open(dir2+'/'+files[idx])
   images2[Nfiles-1-idx,:,:] = hdu2[0].data
   hdu2.close()
   file1 = files[idx].replace('_L3','')
   hdu1 = fits.open(dir1+'/'+file1)
   images1[Nfiles-1-idx,:,:] = hdu1[0].data
   hdu1.close()

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