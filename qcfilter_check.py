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
def qccheck(im, refimages):

    level1 = 0.7

    # first just try something similiar - compare median brightness

    nonzero = np.ma.masked_equal(im,0.0,copy=False)
    med = np.ma.median(nonzero)

    refmeds = []
    for ref in refimages:
        nz = np.ma.masked_equal(ref,0.0,copy=False)
        refmeds.append(np.ma.median(nz))
    refmed = np.array(refmeds).mean()
    if refmed > 0.0:
        rat = med / refmed
    else:
        rat = 1.0

    if rat < level1:
        return 1
    else:
        return 0

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
for file in os.listdir(dir):
    fpath = dir+'/'+file
    try:
        assert(os.path.isfile(fpath))
        hdu = fits.open(fpath)[0]
        images.append(hdu.data)
        files.append(file)
    except Exception as e:
        print(red+f"{e}\nSkipping {file}"+cend)
        pass


#------------------------------------------------------------------------------
# second pass: qc filter

Nimages = len(images)
nx = images[0].shape[0]
ny = images[0].shape[1]

qcflag = np.zeros(Nimages,dtype=int)
for idx in np.arange(Nimages):
    i1 = np.max([0, idx-2])
    i2 = np.min([idx+3, Nimages])
    Nref = i2 - i1
    refimages = np.empty((nx, ny, Nref), dtype = 'float32')
    for ridx in np.arange(Nref):
        refimages[:,:,ridx] = images[ridx+i1]
    qcflag[idx] = qccheck(images[idx], refimages)
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

