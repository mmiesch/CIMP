"""
Driver for the Movie class
"""

from CIMP import Animate as an
from sunpy.net import attrs as a

pcase = 2

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
    title = "LASCO April, 2012"
    outfile = f"/lasco_c3_2012_04_p{pcase}.mp4"
    instrument = a.Instrument.lasco
    detector = a.Detector.c2
    dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/15'
    bgfile = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/background.fts'
    background = 'ratio'
    method = 'none'
    colormap = 'soholasco2'
    scale = (0.01,0.04)

outfile = outdir + '/' + outfile 

x = an.movie(dir, bgfile = bgfile, outfile = outfile, \
             instrument = instrument, detector = detector, \
             cmap = colormap)

x.daymovie(background = background, method = method, scale = scale, \
           title=title)

