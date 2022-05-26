"""
Driver for the Movie class
"""

from CIMP import Animate as an
from sunpy.net import attrs as a

pcase = 10

rmin = 0.0
rmask = None
framedir = None
tolerance = None
diff_ratio = 100.0
resample = None
day = None

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
    scale = (0.5,3.0)

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
    scale = (1.0,2.0)
    resample = 96
    day = '2012-04-15'
    tolerance = 0.02; diff_ratio = 100.0
    #framedir = '/home/mark.miesch/Products/image_processing/frames/debug'

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
    scale = (1.0,2.0)
    resample = 96
    day = '2014-01-17'
    tolerance = 0.01; diff_ratio = 30.0
    framedir = '/home/mark.miesch/Products/image_processing/frames/debug'

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
    scale = (0.0,1.0)

elif pcase == 5:
    title = "Model CME0 pos-30 mgn enhanced"
    outfile = f"/CME0_pos30_p{pcase}_mgn.mp4"
    instrument = 'ModelHAO0'
    detector = 'original'
    dir = '/home/mark.miesch/data/anny/CME0/pos-30/dcmer_030W_bang_0000_fits/tB'
    bgfile = dir+'/frame_0000.fits'
    background = 'ratio'
    method = 'enhance_mgn'
    colormap = 'soholasco2'
    scale = (0.0,1.0)

elif pcase == 6:
    title = "Model CME0 pos-30 fnrgf enhanced"
    outfile = f"/CME0_pos30_p{pcase}_fnrgf.mp4"
    instrument = 'ModelHAO0'
    detector = 'original'
    dir = '/home/mark.miesch/data/anny/CME0/pos-30/dcmer_030W_bang_0000_fits/tB'
    bgfile = dir+'/frame_0000.fits'
    background = 'ratio'
    method = 'enhance_fnrgf'
    colormap = 'soholasco2'
    scale = (0.0,1.0)

elif pcase == 7:
    title = "LASCO April 15, 2012: mgn enhanced"
    outfile = f"/lasco_c3_2012_04_15_p{pcase}_mgn.mp4"
    instrument = a.Instrument.lasco
    detector = a.Detector.c3
    dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/15'
    bgfile = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/background.fts'
    background = 'ratio'
    method = 'enhance_mgn'
    colormap = 'soholasco2'
    #colormap = 'stereocor2'
    rmask = 1.05
    #scale = (0.0,1.0) # v1
    scale = (0.1,1.0)
    tolerance = 0.2; diff_ratio = 30.0
    resample = 96
    day = '2012-04-15'
    framedir = '/home/mark.miesch/Products/image_processing/frames/debug'


elif pcase == 8:
    title = "LASCO April 15, 2012: fnrgf enhanced"
    outfile = f"/lasco_c3_2012_04_15_p{pcase}_fnrgf.mp4"
    instrument = a.Instrument.lasco
    detector = a.Detector.c3
    dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/15'
    bgfile = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/background.fts'
    background = 'ratio'
    method = 'enhance_fnrgf'
    colormap = 'soholasco2'
    #colormap = 'stereocor2'
    rmask = 1.05
    scale = (0.0,1.0)
    tolerance = 0.2; diff_ratio = 30.0
    resample = 96
    day = '2012-04-15'
    framedir = '/home/mark.miesch/Products/image_processing/frames/debug'

elif pcase == 9:
    title = "LASCO Jan 17, 2014: mgn enhanced"
    outfile = f"/lasco_c3_2014_01_17_p{pcase}_mgn.mp4"
    instrument = a.Instrument.lasco
    detector = a.Detector.c3
    dir = '/home/mark.miesch/data/lasco_monthly/c3/2014_01/17'
    bgfile = '/home/mark.miesch/data/lasco_monthly/c3/2014_01/background.fts'
    background = 'ratio'
    method = 'enhance_mgn'
    colormap = 'soholasco2'
    rmask = 1.05
    scale = (0.15,0.9)
    resample = 96
    day = '2014-01-17'
    tolerance = 0.08; diff_ratio = 30.0
    #framedir = '/home/mark.miesch/Products/image_processing/frames/2014_01_17/debug'

elif pcase == 10:
    title = "LASCO Jan 17, 2014: fnrgf enhanced"
    outfile = f"/lasco_c3_2014_01_17_p{pcase}_fnrgf.mp4"
    instrument = a.Instrument.lasco
    detector = a.Detector.c3
    dir = '/home/mark.miesch/data/lasco_monthly/c3/2014_01/17'
    bgfile = '/home/mark.miesch/data/lasco_monthly/c3/2014_01/background.fts'
    background = 'ratio'
    method = 'enhance_fnrgf'
    colormap = 'soholasco2'
    rmask = 1.05
    #scale = (0.0,1.0) v1
    #scale = (0.1,1.0) v2
    scale = (0.2,1.0)
    resample = 96
    day = '2014-01-17'
    tolerance = 0.2; diff_ratio = 30.0
    #framedir = '/home/mark.miesch/Products/image_processing/frames/2014_01_17/debug'

outfile = outdir + '/' + outfile

x = an.movie(dir, bgfile = bgfile, outfile = outfile, \
             instrument = instrument, detector = detector, \
             cmap = colormap)

x.daymovie(background = background, method = method, \
           scale = scale, rmin = rmin, rmax = rmask, title=title, \
           framedir = framedir, tolerance = tolerance, \
           diff_ratio = diff_ratio, resample = resample, day = day)

