"""
Add noise to CME model images
"""

import os

from astropy.io import fits
from skimage.util import random_noise

dir = '/home/mark.miesch/Products/image_processing/ATBD/data/model/CME0_pos30'

indir = dir+'/L2proxy'

#------------------------------------------------------------------------------
ccase = 2
clip = (1.0,4.0)

if ccase == 1:
    outdir = dir + '/L2proxy_gaussian'
else:
    outdir = dir + '/L2proxy_salt'

norm = clip[1] - clip[0]
for file in os.listdir(indir):

    infile = indir+'/'+ file

    hdu = fits.open(infile)
    data = hdu[0].data
    header = hdu[0].header
    hdu.close()
    data = (data - clip[0])/norm
    data = data.clip(min = 0.0, max = 1.0)

    if ccase == 1:
      noisyd = random_noise(data, mode = 'gaussian', var = 0.05)
    else:
      noisyd = random_noise(data, mode = 'salt')

    ndata = noisyd*norm + clip[0]
    outfile = outdir+'/'+file
    hdu_out = fits.PrimaryHDU(ndata, header)
    hdu_out.writeto(outfile, overwrite = True)

#------------------------------------------------------------------------------

