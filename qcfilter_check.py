"""
This script is for development and testing of a QC filter that will flag or remove corrupted images.
"""
red = '\033[91m'
yellow = '\033[93m'
cend = '\033[0m'

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sunpy.visualization.colormaps as cm

from astropy.io import fits

#------------------------------------------------------------------------------
def nzmedian(im):
    nonzero = np.ma.masked_equal(im,0.0,copy=False)
    return np.ma.median(nonzero)

#------------------------------------------------------------------------------
def irange(idx, N, Nref = 5):
    """
    Define a range of indices centered around idx, unless idx is near the edges
    """
    i2 = np.min([idx+3, Nimages])
    if idx < 3:
        range = (0, Nref)
    else:
        range = (i2 - Nref, i2)

    return range

#------------------------------------------------------------------------------
def qc_brightness(med, refmeds):

    qc1 = (0.7,1.3)
    qc2 = (0.5,1.5)

    ref = refmeds.mean()
    if ref > 0.0:
        rat = med/ref
    else:
        rat = 1.0

    if (rat < qc2[0]) | (rat > qc2[1]):
        return 2
    elif (rat < qc1[0]) | (rat > qc1[1]):
        return 1
    else:
        return 0

#------------------------------------------------------------------------------
def qc_diff(images, idx):

    levels = [(1, 0.2, 50), (2, 0.3, 50)]

    refimages = np.array(images)
    ref = np.nanmedian(refimages,axis=0)

    flag = 0
    for lev in reversed(levels):
        d = fits.ImageDataDiff(images[idx], ref, rtol = lev[1])
        if (100*d.diff_ratio > lev[2]):
            flag = lev[0]
            break

    return flag

#------------------------------------------------------------------------------

fig = 2

# default directory for movies
outdir = '/home/mark.miesch/Products/image_processing/movies'

# base directory for frame output
fdir='/home/mark.miesch/Products/image_processing/frames/qc/lasco_2012_04'

framedir = None
skipmovie = False

if fig == 1:

    # This is processed L0.5 data from LASCO/C3

    qclevel = 0
    dir='/home/mark.miesch/data/lasco_monthly/c3/2012_04_processed'
    outfile = 'lasco_proc_qc1.mp4'
    cmap = plt.get_cmap('soholasco2')
    scale = (0.2,1.0)
    framedir=fdir+'/noqc'

elif fig == 2:

    qclevel = 1
    dir='/home/mark.miesch/data/lasco_monthly/c3/2012_04_processed'
    outfile = 'lasco_proc_qc2.mp4'
    cmap = plt.get_cmap('soholasco2')
    scale = (0.2,1.0)
    #framedir=fdir+'/qc2'

else:
    print("pick a valid figure number")
    exit()

#------------------------------------------------------------------------------
# first pass: get images

images = []
files = []
meds = []
for file in os.listdir(dir):
    fpath = dir+'/'+file
    try:
        assert(os.path.isfile(fpath))
        hdu = fits.open(fpath)[0]
        images.append(hdu.data)
        meds.append(nzmedian(hdu.data))
        files.append(file)
    except Exception as e:
        print(red+f"{e}\nSkipping {file}"+cend)
        pass

meds = np.array(meds)

#------------------------------------------------------------------------------
# second pass: apply qc filters

Nimages = len(images)
qcflag = np.zeros(Nimages,dtype=int)
for idx in np.arange(Nimages):
    i = irange(idx,Nimages)
    flag1 = qc_brightness(meds[idx], meds[i[0]:i[1]])
    flag2 = qc_diff(images[i[0]:i[1]], idx - i[0])
    qcflag[idx] = np.max([flag1, flag2])

#------------------------------------------------------------------------------
# print flagged images

for idx in np.arange(Nimages):
    if qcflag[idx] > 0:
        print(f"{idx+1} QC={qcflag[idx]} {files[idx]}")

#------------------------------------------------------------------------------
# remove images that do not pass the qc check

for idx in sorted(np.arange(Nimages), reverse=True):
    if qcflag[idx] >= qclevel:
        del images[idx]
        del files[idx]

#------------------------------------------------------------------------------
# make movie

if skipmovie:
    exit()

fig = plt.figure()
frames = []
for idx in np.arange(len(images)):
    im = plt.figimage(images[idx], cmap=cmap, vmin = scale[0], \
                    vmax = scale[1], origin='lower', resize=True)
    frames.append([im])
    frame = str(len(frames)).zfill(3)

    if framedir is not None:
        plt.savefig(framedir+f"/frame_{frame}.png")

mov = animation.ArtistAnimation(fig, frames, interval = 50, blit = False,
      repeat = True, repeat_delay = 1000)

print(yellow+f"{Nimages} images, {len(frames)} frames"+cend)

mov.save(outdir+'/'+outfile)

