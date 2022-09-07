"""
same as moviemaker but without interpolation onto a regular time grid.  Just take in all files as they are.
"""
#------------------------------------------------------------------------------
import datetime
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import noisegate as ng
import numpy as np
import os
import sunpy.visualization.colormaps as cm
from astropy.io import fits

#------------------------------------------------------------------------------
def brightness_equalization(images, method = None):
    """
    Experiment with brightness equalization for a series of images
    images = 3D numpy array with dimensions (Nt, Nx, Ny)
    """

    if method is None:
        return

    if method == 'mean':
        means = np.nanmean(images, axis=(1,2))
        totmean = np.nanmean(images)
        assert totmean > 0.0, f"ERROR: mean brightness {totmean}"
        for i in np.arange(len(means)):
            images[i,:,:] *= totmean/means[i]

    return

#------------------------------------------------------------------------------
def get_time(header, source):

    if source == 'lasco':
        d = np.array(header['DATE-OBS'].split('/'),dtype='int')
        t = np.rint(np.array(header['TIME-OBS'].split(':'),dtype='float')).astype('int')
        return datetime.datetime(d[0],d[1],d[2],t[0],t[1],t[2])
    else:
        print("ERROR: cannot find time")
        return 0

#------------------------------------------------------------------------------

fig = 6

rootdir = '/home/mark.miesch/Products/image_processing/ATBD'
noisegate = False
framedir = None
beqmethod = None

if fig == 1:
    source = 'lasco'
    title = 'LASCO/C3 April 14-15, 2012'
    cmap = plt.get_cmap('soholasco2')
    dir = rootdir + '/data/lasco_c3/L3_2012_04'
    endfile = 'LASCOC3_L3_2012_04_15_064205.fts'
    duration = 1.0  # duration of movie in days
    scale = (0.0, 1.0)
    outfile = rootdir+'/movies/beta.mp4'

elif fig == 2:
    source = 'lasco'
    title = 'LASCO/C3 April 15-16, 2012'
    cmap = plt.get_cmap('soholasco2')
    dir = rootdir + '/data/lasco_c3/L3_2012_04'
    endfile = 'LASCOC3_L3_2012_04_16_111805.fts'
    duration = 2.0
    scale = (0.0, 1.0)
    outfile = rootdir+'/movies/beta.mp4'

elif fig == 3:
    source = 'lasco'
    title = 'LASCO/C3 Jan 14-16, 2014: no beq'
    cmap = plt.get_cmap('soholasco2')
    dir = rootdir + '/data/lasco_c3/L3_2014_01'
    endfile = 'LASCOC3_L3_2014_01_16_181805.fts'
    duration = 2.0
    scale = (0.0, 0.6)
    beqmethod = None
    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = pdir+f'/movies/lasco2014_beqnorm.mp4'
    framedir = pdir+f'/frames/lasco2014_beqnorm'

elif fig == 4:
    source = 'lasco'
    title = 'LASCO/C3 May 15-17, 2021'
    cmap = plt.get_cmap('soholasco2')
    dir = rootdir + '/data/lasco_c3/L3_2021_05'
    endfile = 'LASCOC3_L3_2021_05_17_013020.fts'
    duration = 2.0
    scale = (0.0, 0.6)
    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/lasco_2021_05_16.mp4'
    framedir = pdir+f'/frames/2021_05_16'

elif fig == 5:
    source = 'lasco'
    title = 'LASCO/C3 May 15-17, 2021'
    cmap = plt.get_cmap('soholasco2')
    dir = rootdir + '/data/lasco_c3/L3_2021_05'
    endfile = 'LASCOC3_L3_2021_05_17_013020.fts'
    duration = 2.0
    scale = (0.0, 0.6)
    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/lasco_2021_05_16_ng.mp4'
    framedir = pdir+f'/frames/2021_05_16_ng'
    noisegate = True

elif fig == 6:
    # try to replicate the atbd
    source = 'lasco'
    title = 'ATBD'
    cmap = plt.get_cmap('soholasco2')
    dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04_processed'
    endfile = 'image167.fts'
    duration = 2.0
    scale = (0.2, 1.0)
    pdir = '/home/mark.miesch/Products/image_processing'
    outfile = rootdir+'/movies/ngcheck.mp4'
    framedir = pdir+f'/frames/ngcheck'
    noisegate = True

else:
    print("pick a valid figure number")
    exit()
#------------------------------------------------------------------------------

# Compile list of valid files in time range of interest,
# with time stamps
dirlist = os.listdir(dir)
flist = list(sorted(dirlist, reverse=True))
idx = flist.index(endfile)

dtmax = datetime.timedelta(days=duration)
dt = datetime.timedelta(days=0.0)

files = []
times = []
while (dt <= dtmax) and (idx < len(flist)):
    hdu = fits.open(dir+'/'+flist[idx])
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

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD debugging
files = list(reversed(os.listdir(dir)))
#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD debugging

Nfiles = len(files)

hdu = fits.open(dir+'/'+files[0])
nx, ny = hdu[0].data.shape
hdu.close()

images = np.zeros((Nfiles,nx,ny), dtype = 'float')
for idx in np.arange(Nfiles):
   hdu = fits.open(dir+'/'+files[idx])
   images[Nfiles-1-idx,:,:] = hdu[0].data
   hdu.close()

#brightness_equalization(images, method = beqmethod)

print(f"MSM {len(images)} {images.shape}")

if noisegate:
    dcube = images.copy()
    images = ng.noise_gate_batch(dcube, cubesize=(18,18,18), model='hybrid',
                                 factor = 4.0, dkfactor = 1.5)

#------------------------------------------------------------------------------
# make movie
fig = plt.figure()
frames = []
for idx in np.arange(Nfiles):
    f = plt.figimage(images[idx,:,:], cmap = cmap, vmin = scale[0], \
            vmax = scale[1], origin='lower', resize=True)
    if title is not None:
        plt.title(title)
    frames.append([f])
    if framedir is not None:
        frame = str(len(frames)).zfill(3)
        plt.savefig(framedir+f"/frame_{frame}.png")

mov = animation.ArtistAnimation(fig, frames, interval = 50, blit = False,
        repeat = True, repeat_delay = 1000)

print(f"Number of valid files = {Nfiles}")
print(f"Number of frames = {len(frames)}")
mov.save(outfile)

