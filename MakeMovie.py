"""
Driver for the Movie class
"""

from CIMP import Animate as an
from sunpy.net import attrs as a

pcase = 16

rmin = 0.0
rmax = None
framedir = None
tolerance = None
diff_ratio = 100.0
resample = None
day = None
morph = False
noisegate = False
downsample = False
clip = None

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
    scale = (0.0,0.1)
    #resample = 96
    #day = '2012-04-15'
    #tolerance = 0.02; diff_ratio = 100.0
    framedir = '/home/mark.miesch/Products/image_processing/frames/debug'

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
    rmax = 1.05
    scale = (0.3,1.0)
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
    rmax = 1.05
    scale = (0.2,1.0)
    #tolerance = 0.2; diff_ratio = 30.0
    tolerance = 0.4; diff_ratio = 30.0
    resample = 96
    day = '2012-04-15'
    #framedir = '/home/mark.miesch/Products/image_processing/frames/debug'

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
    rmax = 1.05
    scale = (0.2,0.9)
    resample = 96
    day = '2014-01-17'
    tolerance = 0.2; diff_ratio = 30.0
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
    rmax = 1.05
    #scale = (0.0,1.0) v1
    #scale = (0.1,1.0) v2
    scale = (0.2,1.0)
    resample = 96
    day = '2014-01-17'
    tolerance = 0.2; diff_ratio = 30.0
    #framedir = '/home/mark.miesch/Products/image_processing/frames/2014_01_17/debug'

elif pcase == 11:
    title = "LASCO April 15, 2012: mgn enhanced+morph pf"
    outfile = f"/lasco_c3_2012_04_15_p{pcase}_mgn_morph.mp4"
    instrument = a.Instrument.lasco
    detector = a.Detector.c3
    dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/15'
    bgfile = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/background.fts'
    background = 'ratio'
    method = 'enhance_mgn'
    colormap = 'soholasco2'
    rmax = 1.05
    scale = (0.15,1.0)
    tolerance = 0.4; diff_ratio = 30.0
    #tolerance = 0.2; diff_ratio = 30.0
    #resample = 96
    #day = '2012-04-15'
    framedir = '/home/mark.miesch/Products/image_processing/frames/2012_04_15/morph'

elif pcase == 12:
    title = "LASCO Jan 17, 2014: mgn enhanced + morph"
    outfile = f"/lasco_c3_2014_01_17_p{pcase}_mgn_morph.mp4"
    instrument = a.Instrument.lasco
    detector = a.Detector.c3
    dir = '/home/mark.miesch/data/lasco_monthly/c3/2014_01/17'
    bgfile = '/home/mark.miesch/data/lasco_monthly/c3/2014_01/background.fts'
    background = 'ratio'
    method = 'enhance_mgn'
    colormap = 'soholasco2'
    rmax = 1.05
    scale = (0.05,0.9)
    resample = 96
    day = '2014-01-17'
    tolerance = 0.2; diff_ratio = 30.0
    framedir = '/home/mark.miesch/Products/image_processing/frames/debug'

elif pcase == 13:
    title = "LASCO April 15, 2012 w/ noisegate"
    outfile = f"/lasco_c3_2012_04_15_p{pcase}_ng.mp4"
    instrument = a.Instrument.lasco
    detector = a.Detector.c3
    dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/15'
    bgfile = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/background.fts'
    background = 'ratio'
    method = 'none'
    #colormap = 'soholasco2'
    colormap = 'stereocor2'
    scale = (0.0,0.1)
    resample = 96
    day = '2012-04-15'
    tolerance = 0.6; diff_ratio = 30.0
    noisegate = True
    framedir = '/home/mark.miesch/Products/image_processing/frames/debug'

elif pcase == 14:
    # stereo L1 data
    title = "STEREO-A Sept 20, 2012"
    #outfile = f"/stereo_a_2012_09_20_p{pcase}_rs72.mp4"
    outfile = f"/stereo_a_2012_09_20_p{pcase}.mp4"
    instrument = a.Instrument.secchi
    detector = a.Detector.cor2
    dir = '/home/mark.miesch/sunpy/data/secchi_cor2/L1/2012/09/20'
    bgfile = '/home/mark.miesch/sunpy/data/secchi_cor2/L1/2012/09/background.fts'
    background = 'ratio'
    method = 'none'
    colormap = 'stereocor2'
    scale = (1.0,1.2)
    rmin = 0.15
    rmax = 1.0
    #resample = 72
    #day = '2012-09-20'
    #tolerance = 0.6; diff_ratio = 30.0
    framedir = '/home/mark.miesch/Products/image_processing/frames/debug'

elif pcase == 15:
    # same as 14 but downsample to 1k x 1k
    title = "STEREO-A Sept 20, 2012"
    outfile = f"/stereo_a_2012_09_20_p{pcase}.mp4"
    instrument = a.Instrument.secchi
    detector = a.Detector.cor2
    dir = '/home/mark.miesch/sunpy/data/secchi_cor2/L1/2012/09/20'
    bgfile = '/home/mark.miesch/sunpy/data/secchi_cor2/L1/2012/09/background.fts'
    background = 'ratio'
    method = 'none'
    colormap = 'stereocor2'
    scale = (1.0,1.3)
    rmin = 0.16
    rmax = 1.0
    resample = 72
    day = '2012-09-20'
    downsample = True
    framedir = '/home/mark.miesch/Products/image_processing/frames/2012_09_20/noenhance'

elif pcase == 16:
    # Try an mgn enhanced version of the stereo L1 data
    title = "STEREO-A Sept 20, 2012"
    outfile = f"/stereo_a_2012_09_20_p{pcase}_mgn.mp4"
    instrument = a.Instrument.secchi
    detector = a.Detector.cor2
    dir = '/home/mark.miesch/sunpy/data/secchi_cor2/L1/2012/09/20'
    bgfile = '/home/mark.miesch/sunpy/data/secchi_cor2/L1/2012/09/background.fts'
    background = 'ratio'
    method = 'enhance_mgn'
    colormap = 'soholasco2'
    clip = (1, 1.3) # clipping the ratio image
    scale = (0.1,1.0)
    rmin = 0.16
    rmax = 1.0
    resample = 72
    day = '2012-09-20'
    downsample = True
    framedir = '/home/mark.miesch/Products/image_processing/frames/2012_09_20/mgn'

elif pcase == 17:
    # Same as 16 but no initial point filter 
    title = "STEREO-A Sept 20, 2012"
    outfile = f"/stereo_a_2012_09_20_p{pcase}_mgn_nopf.mp4"
    instrument = a.Instrument.secchi
    detector = a.Detector.cor2
    dir = '/home/mark.miesch/sunpy/data/secchi_cor2/L1/2012/09/20'
    bgfile = '/home/mark.miesch/sunpy/data/secchi_cor2/L1/2012/09/background.fts'
    background = 'ratio'
    method = 'enhance_mgn'
    #colormap = 'stereocor2'
    colormap = 'soholasco2'
    clip = (1, 1.3) # clipping the ratio image
    scale = (0.1,1.0)
    rmin = 0.15
    rmax = 1.0
    resample = 72
    day = '2012-09-20'
    downsample = True
    framedir = '/home/mark.miesch/Products/image_processing/frames/2012_09_20/mgn_nopf'


outfile = outdir + '/' + outfile

x = an.movie(dir, bgfile = bgfile, outfile = outfile, \
             instrument = instrument, detector = detector, \
             cmap = colormap)

x.daymovie(background = background, method = method, clip = clip, \
           scale = scale, rmin = rmin, rmax = rmax, title=title, \
           framedir = framedir, tolerance = tolerance, \
           diff_ratio = diff_ratio, resample = resample, day = day, \
           noisegate = noisegate, downsample = downsample)

