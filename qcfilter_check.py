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

fig = 1

# default directory for movies
outdir = '/home/mark.miesch/Products/image_processing/movies'

# base directory for frame output
fdir='/home/mark.miesch/Products/image_processing/frames/qc/lasco_2012_04'


if fig == 1:

    # This is processed L0.5 data from LASCO/C3

    dir='/home/mark.miesch/data/lasco_monthly/c3/2012_04_processed'
    outfile = 'lasco_proc_qc1.mp4'
    cmap = plt.get_cmap('soholasco2')
    scale = (0.2,1.0)
    framedir=fdir+'/noqc'

else:
    print("pick a valid figure number")
    exit()

#------------------------------------------------------------------------------

images = []
for file in os.listdir(dir):
    fpath = dir+'/'+file
    try:
        assert(os.path.isfile(fpath))
        hdu = fits.open(fpath)[0]
        images.append(hdu.data)
    except Exception as e:
        print(red+f"{e}\nSkipping {file}"+cend)
        pass

Nimages = len(images)

#------------------------------------------------------------------------------
# make movie
fig = plt.figure()
frames = []
for idx in np.arange(Nimages):
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

