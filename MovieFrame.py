"""
Generate a single movie frame, using the same procedure as in MakeMovie
Intended for debugging
"""

import matplotlib.pyplot as plt

from CIMP import Animate as an
from CIMP import Snapshot as snap
from sunpy.net import attrs as a

pcase = 4

# default directory for movies
outdir = '/home/mark.miesch/Products/image_processing/movies'

if pcase == 1:
    title = "Model CME0 pos-30"
    outfile = f"/CME0_pos30_p{pcase}.mp4"
    instrument = 'ModelHAO0'
    detector = 'original'
    dir = '/home/mark.miesch/data/anny/CME0/pos-30/dcmer_030W_bang_0000_fits/tB'
    bgfile = dir+'/frame_0000.fits'
    background = 'ratio'
    method = 'none'
    colormap = 'soholasco2'
    scale = (0.0,1.0)

elif pcase == 2:
    title = "LASCO April 15, 2012"
    outfile = f"/lasco_c3_2012_04_15_p{pcase}.mp4"
    instrument = a.Instrument.lasco
    detector = a.Detector.c3
    dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/15'
    bgfile = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/background.fts'
    background = 'ratio'
    method = 'none'
    #colormap = 'soholasco2'
    colormap = 'stereocor2'
    scale = (0.0,1.0)

elif pcase == 3:
    title = "LASCO Jan 17, 2014"
    outfile = f"/lasco_c3_2014_01_17_p{pcase}.mp4"
    instrument = a.Instrument.lasco
    detector = a.Detector.c3
    dir = '/home/mark.miesch/data/lasco_monthly/c3/2014_01/17'
    bgfile = '/home/mark.miesch/data/lasco_monthly/c3/2014_01/background.fts'
    background = 'ratio'
    method = 'none'
    colormap = 'soholasco2'
    #colormap = 'stereocor2'
    #scale = (1.,2)
    scale = (0.,1.)

    file1 = '33385478.fts'
    file2 = '33385479.fts'

    #file1 = '33385479.fts'
    #file2 = '33385480.fts'

# subset of simulation data for testing & debugging
elif pcase == 4:
    title = "Testing"
    outfile = f"/testing.mp4"
    instrument = 'ModelHAO0'
    detector = 'original'
    dir = '/home/mark.miesch/data/anny/testing'
    bgfile = dir+'/frame_0000.fits'
    background = 'ratio'
    method = 'enhance mgn'
    colormap = 'soholasco2'
    cliprange = 'image'
    scale = (0.0,1.0)

    file1 = 'frame_0008.fits'
    file2 = 'frame_0009.fits'

outfile = outdir + '/' + outfile 

m = an.movie(dir, bgfile = bgfile, outfile = outfile, \
             instrument = instrument, detector = detector, \
             cmap = colormap)

#-----------------------------------------------------------------


fpath = dir+'/'+file1
x1 = snap.snapshot(file = fpath, bgfile = bgfile, \
    instrument = instrument, detector = detector)
x1.background_ratio(rescale = False)
x1.enhance(detail='mgn')
map1 = x1.map()

fpath = dir+'/'+file2
x2 = snap.snapshot(file = fpath, bgfile = bgfile, \
    instrument = instrument, detector = detector)
x2.background_ratio(rescale = False)
x2.enhance(detail='mgn')
map2 = x2.map()

#-----------------------------------------------------------------

fig = plt.figure(figsize=[16,8])

ax = fig.add_subplot(1,2,1,projection=map1)
map1.plot(cmap = colormap, vmin = scale[0], vmax = scale[1], \
          autoalign = True)

ax = fig.add_subplot(1,2,2,projection=map2)
map2.plot(cmap = colormap, vmin = scale[0], vmax = scale[1], \
          autoalign = True)

plt.show()