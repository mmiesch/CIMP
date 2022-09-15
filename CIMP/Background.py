"""
CIMP Background module
"""

import numpy as np
import os
import sys

from astropy.io import fits

# for warning / error statements; print red, yellow text to terminal
red = '\033[91m'
yellow = '\033[93m'
cend = '\033[0m'

class background:
    """
    The background class is used to compute a background from a series of     coronagraph images.  Typically the series will span 30 days of images     but it will work for other intervals.  The general procedure followed     for computing a background image is:
    1. Compute a daily median image for all images available for a given     day
    2. From the daily median images, construct a background image by     select the minimum positive (non-zero) value for each pixel.
    """

    def __init__(self, dir):
        """
        dir is the directory containing the monthly data.  For now it is assumed that this directory contains one subdirectory for each day of images.  So, for a month of data, there would be 30 subdirectories, each containing all the coronagraph images for a given day.
        """

        self.dir = dir

    def daily_medians(self, normalize = False, lasco_correction = False):
        """
        This loops through the subdirectories for each day and computes a daily median image.  The image is saved in the same directory with the name `daily_median.fits`
        """
        for subdir in os.listdir(self.dir):
            if os.path.isdir(self.dir+'/'+subdir):
                self.compute_median(subdir, normalize, lasco_correction)

    def compute_median(self, daydir, normalize, lasco_correction):
        """
        Given the name of a directory that contains images for a given day, daydir, this method will loop over all the files in that directory to compute a median image.
        """
        Nimages = 0

        ddir = self.dir + '/' + daydir + '/'

        # lasco exposure time correction
        sys.path.append('/home/mark.miesch/Products/image_processing/NRL_python')
        from lasco_utils import get_exp_factor

        caldir = '/home/mark.miesch/data/lasco_cal/idl/expfac/data'

        # for now use the header of the first valid file as a basis for 
        # creating the header for the daily median image file.  So, it will 
        # have the date of the first image - we may want to change this in 
        # the future
        header0 = None

        nx = None
        ny = None

        # first pass: count the number of valid images to process
        # make sure they all hve the same resolution
        files = []
        for file in os.listdir(ddir):
            try:
                assert("median" not in file)
                hdu = fits.open(ddir+file)
                if nx is None:
                    nx = hdu[0].header['NAXIS1']
                else:
                    assert(nx == hdu[0].header['NAXIS1'])

                if ny is None:
                    ny = hdu[0].header['NAXIS1']
                else:
                    assert(ny == hdu[0].header['NAXIS1'])

                etime = hdu[0].header['EXPTIME']

                if header0 is None:
                    header0 = hdu[0].header

                files.append(file)

                hdu.close()

            except:
                print(yellow+f"skipping {file}"+cend)

        Nimages = len(files)
        print(f"{daydir} : {Nimages} : {nx} {ny}")

        assert(Nimages > 0), \
            red+f"ERROR: No valid images for {daydir}"+cend

        # second pass: fill numpy array
        x = np.empty((nx,ny,Nimages),dtype='float32')
        x[:,:,:] = np.nan

        idx = 0
        for file in files:
            hdu= fits.open(ddir+file)
            if lasco_correction:
                etime = hdu[0].header['EXPTIME']
                fac, bias = get_exp_factor(hdu[0].header, dir0 = caldir)
                x[:,:,idx] = (hdu[0].data - bias) / (fac*etime)
                print(yellow+f"Normalized, corrected {etime} {fac} {bias}"+cend)
            elif normalize:
                etime = hdu[0].header['EXPTIME']
                x[:,:,idx] = hdu[0].data / etime
                print(yellow+f"Normalized by exposure time {etime}"+cend)
            else:
                x[:,:,idx] = hdu[0].data

            idx += 1
            hdu.close()

        # compute median image
        med_im = np.nanmedian(x, axis=2)

        # record the number of images used in the header
        header0['NIMAGES'] = Nimages

        # This was causing problems with LASCO with regard to unprintable characters
        del header0['HISTORY']

        hdu_out = fits.PrimaryHDU(data = med_im, header = header0)
        hdu_out.writeto(ddir+'daily_median.fts', overwrite = True)

    def minimize_medians(self):
        """
        This is intended to be executed after daily_medians.
        It takes the daily median file from all the subdirectories
        and computes a background based on the minimum non-zero value
        in each pixel.
        """

        # as before, choose the header in the first file for now
        self.header = None
        self.background = None

        for subdir in os.listdir(self.dir):
            ddir = self.dir+'/'+subdir
            if os.path.isdir(ddir):
                file = ddir+'/'+ 'daily_median.fts'
                try:
                    hdu = fits.open(file)[0]
                    if self.background is None:
                        self.header = hdu.header
                        self.background = hdu.data
                        self.zero = self.background <= 0.0
                    else:
                        np.fmin(self.background, hdu.data, \
                            where = hdu.data > 0.0,
                            out = self.background)
                        self.background = np.where(self.background <= 0.0, \
                                          hdu.data, self.background)
                except:
                    print(yellow+f"skipping {subdir}"+cend)

        mask = self.background < 0.0
        self.background[mask] = 0.0

    def write_background(self):
        bg_file = self.dir + '/' + 'background.fts'
        hdu_out = fits.PrimaryHDU(self.background, self.header)
        hdu_out.writeto(bg_file, overwrite = True)

