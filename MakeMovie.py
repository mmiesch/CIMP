"""
Driver for the Movie class
"""

from CIMP import Animate as an

pcase = 1

outdir = '/home/mark.miesch/Products/image_processing/movies'

if pcase == 1:
    instrument = 'ModelHAO0'
    detector = 'original'
    dir = '/home/mark.miesch/data/anny/CME0/pos-30/dcmer_030W_bang_0000_fits/tB'
    bgfile = dir+'/frame_0000.fits'
    colormap = 'soholasco2'

outfile = outdir + f"/cimp_p{pcase}.mp4"

x = an.movie(dir, bgfile = bgfile, outfile = outfile, \
             instrument = instrument, detector = detector, cmap = colormap)

x.daymovie()

