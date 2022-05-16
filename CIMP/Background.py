"""
CIMP Background module
"""

import numpy as np
import os

from sunpy.io import fits

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

        # for now use the header of the first valid file as a basis for 
        # creating the header for the daily median image file.  So, it will 
        # have the date of the first image - we may want to change this in 
        # the future
        header0 = None

        nx = None
        ny = None

        # first pass: count the number of valid images to process
        # make sure they all hve the same resolution
        for file in os.listdir(ddir):
            try:
                assert("median" not in file)
                data, header = fits.read(ddir+file)[0]
                if nx is None:
                    nx = header['NAXIS1']
                else:
                    assert(nx == header['NAXIS1'])

                if ny is None:
                    ny = header['NAXIS1']
                else:
                    assert(ny == header['NAXIS1'])

                if header0 is None:
                    header0 = header

                Nimages += 1

            except:
                print(yellow+f"skipping {file}"+cend)

        print(f"{daydir} : {Nimages} : {nx} {ny}")

        assert(Nimages > 0), \
            red+f"ERROR: No valid images for {daydir}"+cend

        # second pass: fill numpy array
        x = np.empty((nx,ny,Nimages),dtype='float32')
        x[:,:,:] = np.nan

        idx = 0
        for file in os.listdir(ddir):
            try:
                assert("median" not in file)
                data, header = fits.read(ddir+file)[0]
                assert(nx == header['NAXIS1'])
                assert(ny == header['NAXIS2'])
                x[:,:,idx] = data
                idx += 1
            except:
                pass

        # compute median image
        med_im = np.nanmedian(x, axis=2)

        # record the number of images used in the header
        header0['NIMAGES'] = Nimages

        fits.write(ddir+'daily_median.fts', med_im, header0, overwrite = True)

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
                    data, header = fits.read(file)[0]
                    if self.background is None:
                        self.header = header
                        self.background = data
                        self.zero = self.background <= 0.0
                    else:
                        np.fmin(self.background, data, \
                            where = data > 0.0,
                            out = self.background)
                        self.background = np.where(self.background <= 0.0, \
                                          data, self.background)
                except:
                    print(yellow+f"skipping {subdir}"+cend)

        mask = self.background < 0.0
        self.background[mask] = 0.0

    def write_background(self):
        bg_file = self.dir + '/' + 'background.fts'
        fits.write(bg_file, self.background, self.header, overwrite = True)

