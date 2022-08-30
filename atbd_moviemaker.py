"""
Make a movie by reading in L3 data files, removing corrupted images, applying a noise-gate filter, and nearest-neighbor interpolation on to a regular time grid.
"""
#------------------------------------------------------------------------------
import datetime
import matplotlib.pyplot as plt
import noisegate as ng
import numpy as np
import os
from astropy.io import fits

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

fig = 1

rootdir = '/home/mark.miesch/Products/image_processing/ATBD'

if fig == 1:
    source = 'lasco'
    dir = rootdir + '/data/lasco_c3/L3_2012_04'
    endfile = 'LASCOC3_L3_2012_04_15_064205.fts'
    duration = 1.0  # duration of movie in days
    Nframes = 6 # number of movie frames
    scale = (0.2, 1.0)
    outfile = rootdir+'/movies/beta.mp4'

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
# define regular time grid for movie and identify associated files
# effectively this is nearest-neighbor interpolation, but more 
# memory efficient than reading in all images.

Nfiles = len(files)
tin = np.zeros(Nfiles, dtype='float')
for idx in np.arange(Nfiles):
    tin[idx] = (times[0] - times[idx]).total_seconds()

tgrid = np.linspace(0.0, tin.max(), num = Nframes)

print(tin)
fgrid = []
for t in tgrid:
    diff = np.abs(t-tin)
    idx = np.where(diff == diff.min())[0][0]
    fgrid.append(files[idx])
    print(f"{t} {idx} {tin[idx]} {files[idx]}")

#------------------------------------------------------------------------------
# make movie
#fig = plt.figure()
#frames = []
#for idx in np.arange(Nimages):
#    im = plt.figimage(images[idx], cmap=self.cmap, vmin = scale[0], \
#                    vmax = scale[1], origin='lower', resize=True)
#    if title is not None:
#        plt.title(title)
#    frames.append([im])
#    frame = str(len(frames)).zfill(3)
#    if framedir is not None:
#        plt.savefig(framedir+f"/frame_{frame}.png")
#    mov = animation.ArtistAnimation(fig, frames, interval = 50, blit = False,
#          repeat = True, repeat_delay = 1000)
#    print(f"Number of frames = {len(frames)}")
#    mov.save(outfile)

