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

#------------------------------------------------------------------------------
def nzmedian(im):
    """
    median of nonzero pixels
    """
    nonzero = np.ma.masked_equal(im,0.0,copy=False)
    return np.ma.median(nonzero)

#------------------------------------------------------------------------------
def qc_brightness(images, idx0 = 0):
    """
    QC filter based on changes in median image brightness
    """

    qc1 = (0.7,1.3)
    qc2 = (0.5,1.5)

    N = len(images)
    refmeds = np.zeros(N)
    for idx in np.arange(N):
        refmeds[idx] = nzmedian(images[idx])
    ref = np.nanmedian(refmeds)

    if ref > 0.0:
        rat = refmeds[idx0]/ref
    else:
        rat = 1.0

    if (rat < qc2[0]) | (rat > qc2[1]):
        print(f"qc_brightness flag 2 {rat}")
        return 2
    elif (rat < qc1[0]) | (rat > qc1[1]):
        print(f"qc_brightness flag 1 {rat}")
        return 1
    else:
        return 0

#------------------------------------------------------------------------------
def qc_diff(images, idx0 = 0):
    """
    QC filter based on direct image comparisons
    """

    levels = [(1, 0.2, 50), (2, 0.3, 50)]

    refimages = np.array(images)
    ref = np.nanmedian(refimages,axis=0)

    flag = 0
    for lev in reversed(levels):
        d = fits.ImageDataDiff(images[idx0], ref, rtol = lev[1])
        if (100*d.diff_ratio > lev[2]):
            print(f"qc_diff flag {lev[0]} {lev[1]} {100*d.diff_ratio}")
            flag = lev[0]
            break

    return flag

#------------------------------------------------------------------------------

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
        self.outdir = outdir

        # generate output filename
        filename = os.path.basename(infile).split('_')
        filename.insert(1,'L3')
        self.outfile = outdir+'/'+'_'.join(filename)

    def process(self, rmin = 0.0, rmax = np.inf, clip = None):
        # Proposed L3 pipeline
        # not including noise reduction and QC

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

    def qcfilter(self, Nref = 5):
        """
        Apply QC filter based on Nref-1 previous reference images.
        Assume for now files are ordered alphabetically via time stamp.
        """

        # get existing files in L3 directory
        ofile = os.path.basename(self.outfile)
        dirlist = os.listdir(self.outdir)
        dirlist.append(ofile)
        slist = list(sorted(dirlist, reverse=True))
        idx = slist.index(ofile) + 1

        images = [self.data]
        rfiles = []
        while (len(images) < Nref) and (idx < len(slist)):
            fpath = self.outdir+'/'+slist[idx]
            hdu = fits.open(fpath)
            try:
                flag = hdu[0].header['L3QCFLAG']
            except:
                flag = 0
            if flag < 2:
                images.append(hdu[0].data)
                rfiles.append(slist[idx])
            idx += 1

        if len(images) < 2:
            self.header['L3QCFLAG'] = 0
            return

        flag1 = qc_brightness(images)
        flag2 = qc_diff(images)
        self.header['L3QCFLAG']= np.max([flag1, flag2])
        return

    def write(self):

        assert(self.outfile != self.infile)
        hdu_out = fits.PrimaryHDU(self.data, self.header)
        hdulist = fits.HDUList([hdu_out])
        hdulist.writeto(self.outfile, overwrite = True)
        hdulist.close()