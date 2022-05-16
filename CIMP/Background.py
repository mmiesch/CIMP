"""
CIMP Background module
"""

import numpy as np
import os
import sunpy.io

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

    def daily_medians(self):
        """
        This loops through the subdirectories for each day and computes a daily median image.  The image is saved in the same directory with the name `daily_median.fits`
        """

        for subdir in os.listdir(self.dir):
            if os.path.isdir(self.dir+'/'+subdir):
                self.compute_median(subdir)


    def compute_median(self, daydir):
        """
        Given the name of a directory that contains images for a given day, daydir, this method will loop over all the files in that directory to compute a median image.
        """
        Nimages = 0

        ddir = self.dir + '/' + daydir + '/'

        nx = None
        ny = None

        # first pass: count the number of valid images to process
        # make sure they all hve the same resolution
        for file in os.listdir(ddir):
            try:
                data, header = sunpy.io.fits.read(ddir+file)[0]
                if nx is None:
                    nx = header['NAXIS1']
                else:
                    assert(nx == header['NAXIS1'])

                if ny is None:
                    ny = header['NAXIS1']
                else:
                    assert(ny == header['NAXIS1'])

                Nimages += 1

            except:
                pass

        print(f"{daydir} : {Nimages} : {nx} {ny}")

        assert(Nimages > 0), \
            red+f"ERROR: No valid images for {daydir}"+cend

        # second pass: fill numpy array
        x = np.empty((nx,ny,Nimages),dtype='float32')
        x[:,:,:] = np.nan

        idx = 1
        for file in os.listdir(ddir):
            try:
                data, header = sunpy.io.fits.read(ddir+file)[0]
                assert(nx == header['NAXIS1'])
                assert(ny == header('NAXIS2'))
            except:
                pass

