"""
This module contains prototype code for producing CCOR-1 and CCOR-2 L3 data products.
"""

import noisegate as ng
import numpy as np
import os

from astropy.io import fits
from CIMP import Enhance

# for warning / error statements; print red, yellow text to terminal
red = '\033[91m'
yellow = '\033[93m'
cend = '\033[0m'

class l3proc_ng:
    """
    Class for CCOR L3 data processing with noise-gate.  The use of noise-gate filtering requires the analysis of N image
    """

    def __init__(self, infile, outdir, Nimages = 18, rmin = 0.0, rmax = np.inf):
        """
        infile: This is intended to represent a new L1b (CCOR-1) or L2 (CCOR-2) input file that has been created as part of a real-time operational pipeline.

        outdir: The output directory where the L3 data should be written

        Nimages: the number of images in the image sequence.  This is used by the noise-gate filter to estimate the noise spectrum.

        rmin, rmax: FOV (normalized for minimum extent of image axes) for mask_annulus

        """

        self.infile = os.path.basename(infile)
        self.Nimages = Nimages

        # parent directory containing L1b (CCOR-1) or L2 (CCOR-1) image files
        self.indir = os.path.dirname(infile)
        self.outdir = outdir
        assert(self.indir != self.outdir)

        # loop over files to obtain the most recent N files that pass the
        # QC filter.  In the actual operational implementation this would look
        # different.  However, for now, take advantange of my standard naming
        # convention for L2 proxy data.  In this convention, the filename is
        # includes the time stamp in such a way that an alphabetic ordering 
        # of the files will also be a chronological ordering.  So, start with
        # the last file in an alphabetical ordering and work backwards until
        # you get to N=18

        flist = list(sorted(os.listdir(self.indir)))

        try:
            idx = flist.index(self.infile)
        except ValueError:
            print(red+"ERROR: infile not found"+cend)
            exit()

        # get resolution from first file
        hdu = fits.open(self.indir+'/'+flist[idx])
        self.Nx, self.Ny = hdu[0].data.shape
        hdu.close()

        # resolution of downsampled image
        self.nx, self.ny = int(self.Nx/2), int(self.Ny/2)

        # mask annulus, downsample, and QC check during the initial read

        self.files = []
        self.headers = []
        self.images = np.zeros((self.Nimages, self.nx, self.ny), dtype='float32')
        while len(self.files) < self.Nimages:
            hdu = fits.open(self.indir+'/'+flist[idx])
            data = hdu[0].data
            header = hdu[0].header
            Enhance.mask_annulus(data, rmin = rmin, rmax = rmax)
            data = Enhance.downsample(data)
            assert(data.shape == (self.nx, self.ny))
            header['NAXIS1'] = self.nx
            header['NAXIS2'] = self.ny
            header['CRPIX1'] /= 2
            header['CRPIX2'] /= 2
            header['CDELT1'] *= 2
            header['CDELT2'] *= 2
            self.images[len(self.files),:,:] = data
            self.files.append(flist[idx])
            self.headers.append(hdu[0].header)
            hdu.close()
            idx -= 1
