
"""
Script to create movies showing images before and after processing
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
    elif source == 'model':
        d = header['DATE_OBS']
        return datetime.datetime.fromisoformat(d)
    else:
        print("ERROR: cannot find time")
        return 0

#------------------------------------------------------------------------------
def noisegate(images, cubesize=(18,18,18), factor = 6.0):

    #nap = int((2*cubesize[0])/3)
    nap = 14
    nt, nx, ny = images.shape
    dcube = np.zeros((nt+nap, nx, ny), dtype = 'float')
    dcube[0:nt,:,:] = images
    dcube[nt:,:,:] = np.flip(images[nt-nap:,:,:], axis = 0)

    print(f"Applying noisegate {images.shape} {dcube.shape}")

    dcubeng = ng.noise_gate_batch(dcube, cubesize=cubesize, model='constant',
                factor = factor)
    return dcubeng[nap:-nap,:,:]

#------------------------------------------------------------------------------
fig = 6

rootdir = '/home/mark.miesch/Products/image_processing/ATBD'
ngflag = True
framedir = None
factor = 6.0

if fig == 1:
    source = 'lasco'
    dir2 = rootdir + '/data/lasco_c3/L3_2014_01'
    cmap2 = plt.get_cmap('soholasco2')
    endfile = 'LASCOC3_L3_2014_01_17_053005.fts'
    duration = 2.0
    #scale2 = (0.2, 0.6)
    scale2 = (0.24, 0.6)

    dir1 = rootdir + '/data/lasco_c3/L2proxy_2014_01'
    cmap1 = plt.get_cmap('stereocor2')
    scale1 = (1.0, 1.3)

    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/lasco_2014_01_16_ba.mp4'
    #framedir = pdir+f'/frames/2014_01_16_ba'
    framedir = pdir+f'/frames/beta'

elif fig == 2:
    source = 'stereo'
    dir2 = rootdir + '/data/stereo_a/L3_2012_09'
    cmap2 = plt.get_cmap('soholasco2')
    endfile = 'STEREOA_L3_2012_09_17_022400.fts'
    duration = 2.0
    scale2 = (0.0, 1.0)

    dir1 = rootdir + '/data/stereo_a/L2proxy_2012_09'
    cmap1 = plt.get_cmap('stereocor2')
    scale1 = (1.0, 1.4)

    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/stereo_2012_09_16_ba.mp4'
    framedir = pdir+f'/frames/stereo_2012_09_16_ba'

elif fig == 3:
    # abandoned attempt to work with lasco L1
    source = 'lascoL1'
    dir2 = rootdir + '/data/lasco_c3/L3_L1_2014_01'
    cmap2 = plt.get_cmap('soholasco2')
    endfile = 'LASCOC3_L3_L1_2014_01_17_050449.fts'
    duration = 2.0
    scale2 = (0.0, 0.6)
    factor = 6.0

    dir1 = rootdir + '/data/lasco_c3/L2proxy1_2014_01'
    cmap1 = plt.get_cmap('stereocor2')
    scale1 = (1.0, 1.2)

    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/lascoL1_2014_01_16_ba.mp4'
    #framedir = pdir+f'/frames/2021_05_16_ba'

elif fig == 4:
    # main result for lasco 2021 data
    source = 'lasco'
    dir2 = rootdir + '/data/lasco_c3/L3_2021_05'
    cmap2 = plt.get_cmap('soholasco2')
    endfile = 'LASCOC3_L3_2021_05_23_223007.fts'
    duration = 2.0
    scale2 = (0.0, 0.6)
    factor = 6.0

    dir1 = rootdir + '/data/lasco_c3/L2proxy_2021_05'
    cmap1 = plt.get_cmap('stereocor2')
    scale1 = (1.0, 1.3)

    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/lasco_2021_05_23_ba.mp4'
    framedir = pdir+f'/frames/2021_05_23_f6'

elif fig == 5:
    # experimenting with lasco 2021 data
    source = 'lasco'
    dir2 = rootdir + '/data/lasco_c3/L3_2021_05'
    cmap2 = plt.get_cmap('soholasco2')
    endfile = 'LASCOC3_L3_2021_05_23_223007.fts'
    duration = 2.0
    scale2 = (0.0, 0.6)
    factor = 10.0
    ngflag = False

    dir1 = rootdir + '/data/lasco_c3/L2proxy_2021_05'
    cmap1 = plt.get_cmap('stereocor2')
    scale1 = (1.0, 1.3)

    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/lasco_2021_05_23_nong.mp4'
    framedir = pdir+f'/frames/2021_05_23_nong'

elif fig == 6:
    # main result for lasco 2021 data
    source = 'lasco'
    dir2 = rootdir + '/data/lasco_c3/L3_2012_04'
    cmap2 = plt.get_cmap('soholasco2')
    endfile = 'LASCOC3_L3_2012_04_16_111805.fts'
    duration = 2.0
    scale2 = (0.0, 0.6)
    factor = 6.0

    dir1 = rootdir + '/data/lasco_c3/L2proxy_2012_04'
    cmap1 = plt.get_cmap('stereocor2')
    scale1 = (1.0, 1.3)

    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/lasco_2012_04_16.mp4'
    #framedir = pdir+f'/frames/2012_04_16'

elif fig == 7:
    source = 'model'
    dir2 = rootdir + '/data/model/CME0_pos30/L3'
    cmap2 = plt.get_cmap('soholasco2')
    endfile = 'Model0_L3_2010_04_17_234614.fts'
    duration = 4.0
    scale2 = (0.0, 0.6)
    #ngflag = False

    dir1 = rootdir + '/data/model/CME0_pos30/L2proxy'
    cmap1 = plt.get_cmap('stereocor2')
    scale1 = (1.0, 6.0)

    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/CME0.mp4'
    framedir = pdir+f'/frames/CME0'

elif fig == 8:
    source = 'model'
    dir2 = rootdir + '/data/model/CME0_pos30/L3_gaussian'
    cmap2 = plt.get_cmap('soholasco2')
    endfile = 'Model0_L3_2010_04_17_234614.fts'
    duration = 4.0
    scale2 = (0.0, 0.6)
    #ngflag = False

    dir1 = rootdir + '/data/model/CME0_pos30/L2proxy_gaussian'
    cmap1 = plt.get_cmap('stereocor2')
    scale1 = (1.0, 6.0)

    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/CME0_gaussian.mp4'
    framedir = pdir+f'/frames/CME0_gaussian'

elif fig == 9:
    source = 'model'
    dir2 = rootdir + '/data/model/CME0_pos30/L3_salt'
    cmap2 = plt.get_cmap('soholasco2')
    endfile = 'Model0_L3_2010_04_17_234614.fts'
    duration = 4.0
    scale2 = (0.0, 0.6)
    #ngflag = False

    dir1 = rootdir + '/data/model/CME0_pos30/L2proxy_salt'
    cmap1 = plt.get_cmap('stereocor2')
    scale1 = (1.0, 6.0)

    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/CME0_salt.mp4'
    framedir = pdir+f'/frames/CME0_salt'

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

if ngflag:
    images2 = noisegate(images2, factor = factor)
    #nt = 12 # for cubesize of 18
    nt = 14
    images1 = images1[nt:,:,:]

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

