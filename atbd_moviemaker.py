"""
Make a movie by reading in L3 data files, removing corrupted images, applying a noise-gate filter, and nearest-neighbor interpolation on to a regular time grid.
"""
#------------------------------------------------------------------------------
import datetime
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

rootdir = '/home/mark.miesch/Products/image_processing/ATBD/data'

if fig == 1:
    source = 'lasco'
    dir = rootdir + '/lasco_c3/L3_2012_04'
    endfile = 'LASCOC3_L3_2012_04_15_064205.fts'
    duration = 1.0  # duration of movie in days

#------------------------------------------------------------------------------

# load images
dirlist = os.listdir(dir)
files = list(sorted(dirlist, reverse=True))
idx = files.index(endfile)

dtmax = datetime.timedelta(days=duration)
dt = datetime.timedelta(days=0.0)

images = []
times = []
while (dt <= dtmax) and (idx < len(files)):
    hdu = fits.open(dir+'/'+files[idx])
    try:
        flag = hdu[0].header['L3QCFLAG']
    except:
        flag = 0
    if flag < 2:
        t = get_time(hdu[0].header, source)
        if len(times) > 0:
            dt = times[0] - t
        print(f"{files[idx]} {dt}")
        images.append(hdu[0].data)
        times.append(t)
    idx += 1
