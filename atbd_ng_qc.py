"""
This script is designed to do some qc for applying the noisegate filter,
created while writing the CCOR ATBD.  The starting point is a directory
of images in fits files that have already been background-subtracted
and processed with, e.g. OMR and MGN.  The goal is to use the dirmovie()
method of the Animate class to go through these images and remove any
that are corrupted.  The result is a directory of clean images that
can be used for noisegate noise modeling.   As a bonus, it will also
create a movie.
"""

from CIMP import Animate as an
from sunpy.net import attrs as a

pcase = 1

fitsdir = None
framedir = None
denylist = None
tolerance = None
diff_ratio = 100.0

# default directory for movies
outdir = '/home/mark.miesch/Products/image_processing/movies'

if pcase == 1:
    title = "LASCO April 14-16, 2012"
    outfile = f"/lasco_2012_04_atbd_qc.mp4"
    dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04_processed'
    colormap = 'soholasco2'
    #colormap = 'stereocor2'
    scale = (0.2,1.0)
    tolerance = 0.2; diff_ratio = 50.0
    fitsdir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04_qc'
    framedir = '/home/mark.miesch/Products/image_processing/frames/debug'
    denylist = [4, 19, 22, 27, 45, 60, 98, 111, 114, 117, 122, 124, 128, 131, \
                140, 141, 147, 150, 168]
else:
    print("pick a valid pcase")
    exit()

outfile = outdir + '/' + outfile

x = an.movie(dir, outfile = outfile, cmap = colormap)

x.dirmovie(scale = scale, title=title, denylist = denylist, \
           tolerance = tolerance, diff_ratio = diff_ratio, \
           framedir = framedir, fitsdir = fitsdir)

