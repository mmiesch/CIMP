
"""
Script to compare different ng methods for the ATBD
"""

import datetime
import matplotlib.animation as animation
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
        t = np.floor(np.array(header['TIME-OBS'].split(':'),dtype='float')).astype('int')
        return datetime.datetime(d[0],d[1],d[2],t[0],t[1],t[2])
    elif source == 'stereo' or source == 'lascoL1':
        d = header['DATE-OBS']
        return datetime.datetime.fromisoformat(d)
    else:
        print("ERROR: cannot find time")
        return 0

#------------------------------------------------------------------------------
def noisegate(images, cubesize=(12,12,12), factor = 4.0, model = 'constant', \
              dkfactor = 4.0):

    nap = int((2*cubesize[0])/3)
    nt, nx, ny = images.shape
    dcube = np.zeros((nt+nap, nx, ny), dtype = 'float')
    dcube[0:nt,:,:] = images
    dcube[nt:,:,:] = np.flip(images[nt-nap:,:,:], axis = 0)

    print(f"Applying noisegate {images.shape} {dcube.shape}")

    if model == 'hybrid':
        dcubeng = ng.noise_gate_batch(dcube, cubesize=cubesize, model=model,
                    factor = factor, dkfactor = dkfactor)
    else:
        dcubeng = ng.noise_gate_batch(dcube, cubesize=cubesize, model=model,
            factor = factor)

    return dcubeng[nap:-nap,:,:]

#------------------------------------------------------------------------------
fig = 6

rootdir = '/home/mark.miesch/Products/image_processing/ATBD'
ngflag1 = True
framedir = None

if fig == 1:
    source = 'lasco'
    dir1 = rootdir + '/data/lasco_c3/L3_2021_05'
    dir2 = dir1
    cmap1 = plt.get_cmap('soholasco2')
    cmap2 = cmap1
    endfile = 'LASCOC3_L3_2021_05_23_223007.fts'
    duration = 2.0

    model1 = 'hybrid'
    cubesize1=(18,18,18)
    factor1 = 4.0
    dkfactor1 = 6.0
    scale1 = (0.0, 0.6)

    model2 = 'constant'
    cubesize2=(18,18,18)
    factor2 = 6.0
    dkfactor2 = 6.0
    scale2 = (0.0, 0.6)

    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/lasco_2021_05_ng_ch.mp4'
    framedir = pdir+f'/frames/2021_05_ng_ch'

elif fig == 2:
    source = 'lasco'
    dir1 = rootdir + '/data/lasco_c3/L3_2021_05'
    dir2 = dir1
    cmap1 = plt.get_cmap('soholasco2')
    cmap2 = cmap1
    endfile = 'LASCOC3_L3_2021_05_23_223007.fts'
    duration = 2.0

    model1 = 'constant'
    cubesize1=(18,18,18)
    factor1 = 4.0
    dkfactor1 = 6.0
    scale1 = (0.0, 0.6)

    model2 = 'constant'
    cubesize2=(18,18,18)
    factor2 = 6.0
    dkfactor2 = 6.0
    scale2 = (0.0, 0.6)

    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/lasco_2021_05_ng_46.mp4'
    framedir = pdir+f'/frames/2021_05_ng_46'

elif fig == 3:
    source = 'lasco'
    dir1 = rootdir + '/data/lasco_c3/L3_2021_05'
    dir2 = dir1
    cmap1 = plt.get_cmap('soholasco2')
    cmap2 = cmap1
    endfile = 'LASCOC3_L3_2021_05_23_223007.fts'
    duration = 2.0

    ngflag1 = False
    model1 = 'constant'
    cubesize1=(18,18,18)
    factor1 = 4.0
    dkfactor1 = 6.0
    scale1 = (0.0, 0.6)

    model2 = 'constant'
    cubesize2=(18,18,18)
    factor2 = 6.0
    dkfactor2 = 6.0
    scale2 = (0.0, 0.6)

    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/lasco_2021_05_ng.mp4'
    framedir = pdir+f'/frames/2021_05_ng'

elif fig == 4:
    source = 'lasco'
    dir1 = rootdir + '/data/lasco_c3/L3_2021_05'
    dir2 = dir1
    cmap1 = plt.get_cmap('soholasco2')
    cmap2 = cmap1
    endfile = 'LASCOC3_L3_2021_05_23_223007.fts'
    duration = 2.0

    model1 = 'constant'
    cubesize1=(12,12,12)
    factor1 = 6.0
    dkfactor1 = 6.0
    scale1 = (0.0, 0.6)

    model2 = 'constant'
    cubesize2=(18,18,18)
    factor2 = 6.0
    dkfactor2 = 6.0
    scale2 = (0.0, 0.6)

    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/lasco_2021_05_ng_c18.mp4'
    framedir = pdir+f'/frames/2021_05_ng_c18'

elif fig == 5:
    source = 'stereo'
    dir1 = rootdir + '/data/stereo_a/L3_2012_09'
    dir2 = dir1
    cmap1 = plt.get_cmap('soholasco2')
    cmap2 = cmap1
    endfile = 'STEREOA_L3_2012_09_17_022400.fts'
    duration = 2.0

    model1 = 'hybrid'
    cubesize1=(18,18,18)
    factor1 = 4.0
    dkfactor1 = 6.0
    scale1 = (0.0, 1.0)

    model2 = 'constant'
    cubesize2=(18,18,18)
    factor2 = 6.0
    dkfactor2 = 6.0
    scale2 = (0.0, 1.0)

    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/stereo_2012_09_ng_ch.mp4'
    framedir = pdir+f'/frames/stereo_2012_09_ng_ch'

elif fig == 6:
    source = 'stereo'
    dir1 = rootdir + '/data/stereo_a/L3_2012_09'
    dir2 = dir1
    cmap1 = plt.get_cmap('soholasco2')
    cmap2 = cmap1
    endfile = 'STEREOA_L3_2012_09_17_022400.fts'
    duration = 2.0

    model1 = 'hybrid'
    cubesize1=(12,12,12)
    factor1 = 4.0
    dkfactor1 = 6.0
    scale1 = (0.0, 1.0)

    model2 = 'constant'
    cubesize2=(18,18,18)
    factor2 = 6.0
    dkfactor2 = 6.0
    scale2 = (0.0, 1.0)

    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/stereo_2012_09_ng_c18.mp4'
    framedir = pdir+f'/frames/stereo_2012_09_ng_c18'

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

hdu = fits.open(dir2+'/'+files[0])
nx1, ny1 = hdu[0].data.shape
nx2, ny2 = nx1, ny1
hdu.close()

images = np.zeros((Nfiles,nx1,ny1), dtype = 'float')
for idx in np.arange(Nfiles):
   hdu2 = fits.open(dir2+'/'+files[idx])
   images[Nfiles-1-idx,:,:] = hdu2[0].data
   hdu2.close()

if ngflag1:
    images1 = noisegate(images, cubesize=cubesize1, factor = factor1, \
                        model=model1, dkfactor=dkfactor1)
else:
    nt = 12
    images1 = images[nt:,:,:]

images2 = noisegate(images, cubesize=cubesize2, factor = factor2, \
                    model=model2, dkfactor=dkfactor2)

if cubesize1[0] == 12 and cubesize2[0] == 18:
    images1 = images1[4:,:,:]

#------------------------------------------------------------------------------

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,6))
fig.tight_layout(pad=1,rect=(0.01,0.01,.99,.99))

frames = []

for idx in np.arange(images2.shape[0]):

    f1 = ax[0].imshow(images1[idx,:,:],cmap=cmap1, vmin=scale1[0], vmax=scale1[1], \
                 origin='lower')
    ax[0].axis('off')

    f2 = ax[1].imshow(images2[idx,:,:],cmap=cmap2, vmin=scale2[0], vmax=scale2[1], \
                 origin='lower')
    ax[1].axis('off')
    frames.append([f1,f2])
    if framedir is not None:
        frame = str(len(frames)).zfill(3)
        plt.savefig(framedir+f"/frame_{frame}.png")

mov = animation.ArtistAnimation(fig, frames, interval = 50, blit = False,
        repeat = True, repeat_delay = 1000)

print(f"Number of valid files = {Nfiles}")
print(f"Number of frames = {len(frames)}")
mov.save(outfile)

