"""
Driver for the Movie class
"""

from CIMP import Animate as an
from sunpy.net import attrs as a

pcase = 5

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
    cliprange = 'image'
    scale = (1.0,2.0)

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
    cliprange = (1.0,2.0)
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
    #colormap = 'soholasco2'
    colormap = 'stereocor2'
    cliprange = (1.0,2.0)
    scale = (1.0,2.0)

# subset of simulation data for testing & debugging
elif pcase == 4:
    title = "Testing"
    outfile = f"/testing.mp4"
    instrument = 'ModelHAO0'
    detector = 'original'
    dir = '/home/mark.miesch/data/anny/testing'
    bgfile = dir+'/frame_0000.fits'
    background = 'ratio'
    method = 'enhance_mgn'
    colormap = 'soholasco2'
    cliprange = 'image'
    scale = (0.0,1.0)

elif pcase == 5:
    title = "Testing"
    outfile = f"/CME0_pos30_p{pcase}_mgn.mp4"
    instrument = 'ModelHAO0'
    detector = 'original'
    dir = '/home/mark.miesch/data/anny/CME0/pos-30/dcmer_030W_bang_0000_fits/tB'
    bgfile = dir+'/frame_0000.fits'
    background = 'ratio'
    method = 'enhance_mgn'
    colormap = 'soholasco2'
    cliprange = 'image'
    scale = (0.0,1.0)

outfile = outdir + '/' + outfile 

x = an.movie(dir, bgfile = bgfile, outfile = outfile, \
             instrument = instrument, detector = detector, \
             cmap = colormap)

x.daymovie(background = background, method = method, \
           scale = scale, title=title)

