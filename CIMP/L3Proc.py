"""
This module contains prototype code for producing CCOR-1 and CCOR-2 L3 data products.
"""

import noisegate as ng
import numpy as np
import os

from astropy.io import fits
from CIMP import Enhance
from sunkit_image.enhance import mgn

# for warning / error statements; print red, yellow text to terminal
red = '\033[91m'
yellow = '\033[93m'
cend = '\033[0m'

class l3proc:
    """
    Class for CCOR L3 data processing with noise-gate.  The use of noise-gate filtering requires the analysis of N image
    """

    def __init__(self, infile, outdir):
        """
        infile: This is intended to represent a new L1b (CCOR-1) or L2 (CCOR-2) input file that has been created as part of a real-time operational pipeline.

        outdir: The output directory where the L3 data should be written

        rmin, rmax: FOV (normalized for minimum extent of image axes) for mask_annulus

        """

        self.infile = infile

        # generate output filename
        filename = os.path.basename(infile).split('_')
        filename.insert(1,'L3')
        self.outfile = outdir+'/'+'_'.join(filename)


    def process(self, rmin = 0.0, rmax = np.inf, clip = None):
        # Basic L3 pipeline

        hdu = fits.open(self.infile)
        indata = hdu[0].data
        self.header = hdu[0].header

        # median downsample
        self.data = Enhance.downsample(indata)
        self.nx, self.ny = self.data.shape
        self.header['NAXIS1'] = self.nx
        self.header['NAXIS2'] = self.ny
        self.header['CRPIX1'] /= 2
        self.header['CRPIX2'] /= 2
        self.header['CDELT1'] *= 2
        self.header['CDELT2'] *= 2

        # mask annulus
        Enhance.mask_annulus(self.data, rmin = rmin, rmax = rmax)

        # OMR point removal
        self.data = Enhance.omr(self.data, rescaleim = False)

        # clip and rescale
        self.data = Enhance.clip(self.data, min = clip[0], max = clip[1], rescale_output = True)

        # MGN feature enhancement
        self.data = mgn(self.data, h = 0.8, gamma = 1.5)

        # mask annulus again
        Enhance.mask_annulus(self.data, rmin = rmin, rmax = rmax)

        hdu.close()

    def write(self):

        assert(self.outfile != self.infile)
        hdu_out = fits.PrimaryHDU(self.data, self.header)
        hdulist = fits.HDUList([hdu_out])
        hdulist.writeto(self.outfile, overwrite = True)
        hdulist.close()