
"""
The purpose of this script is to process a bunch of images for ingestion into the L3proc processing pipeline.  So the output images are fits files that are intended to mimic CCOR L2 or L1b data.  This mainly means that the background has been subtracted and the images have been normalized by exposure time.
"""
import os
import numpy as np

from astropy.io import fits
from CIMP import Snapshot as snap
from sunpy.net import attrs as a

#------------------------------------------------------------------------------

outroot = '/home/mark.miesch/Products/image_processing/ATBD/data/'

rmin = 0.16
rmax = 1.0

fig = 1

if fig == 1:

    # L0.5 LASCO data
    instrument = a.Instrument.lasco
    detector = a.Detector.c3
    dir='/home/mark.miesch/data/lasco_monthly/c3/2012_04'
    bgfile = dir+'/background.fts'

    outdir = outroot+'lasco_c3/L2proxy_2012_04'

else:
    print("pick a valid figure number")
    exit()

if not os.path.exists(outdir):
    os.mkdir(outdir)

#------------------------------------------------------------------------------
# loop over all directories for the month

for d in os.listdir(dir):
    day = dir+'/'+d
    print(day)
    if os.path.isdir(day):
        for file in os.listdir(day):
            fpath = day+'/'+file
            try:
                assert(os.path.isfile(fpath))
                assert("median" not in file)
                assert(fpath != bgfile)
                x = snap.snapshot(file = fpath, bgfile = bgfile, \
                    instrument = instrument, detector = detector, \
                    normalize = False)
                t = x.time.datetime
                tstamp = f"{t.year}_{str(t.month).zfill(2)}_{str(t.day).zfill(2)}_{str(t.hour).zfill(2)}{str(t.minute).zfill(2)}{str(t.second).zfill(2)}"
                outfile=outdir+'/'+f"{x.instrument}{x.detector}_{tstamp}"
                x.mask_background(rmin = rmin, rmax = rmax, nonzero = True)
                x.background_ratio(rescale=False)
                header0 = x.header
                # needed for lasco data because the header includes an
                # unprintable character that confuses astropy
                try:
                    del header0['HISTORY']
                except:
                    pass
                hdu_out = fits.PrimaryHDU(x.data,header0)
                hdu_out.writeto(outfile, overwrite = True)
            except Exception as e:
                print(f"{e}\nSkipping {file}")
                pass




#------------------------------------------------------------------------------
