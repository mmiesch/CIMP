"""
CIMP Event module
"""

import astropy.io
import astropy.units as u
import datetime
import logging
from astropy.utils import data_info
import numpy as np
import os
import sunpy.map

from astropy.io import fits
from CIMP import Enhance
from io import RawIOBase
from sunpy.net import Fido
from sunpy.net import attrs as a
from skimage import exposure

# for warning / error statements; print red, yellow text to terminal
red = '\033[91m'
yellow = '\033[93m'
cend = '\033[0m'

testsnap = {
    1: {
    'instrument': a.Instrument.lasco,
    'detector': a.Detector.c3,
    'dir': '/home/mark.miesch/data/lasco_monthly/c3/2012_04/',
    'file': '/15/32296650.fts',
    'bgfile': 'background.fts'
    },
    2: {
    'instrument': a.Instrument.lasco,
    'detector': a.Detector.c3,
    'dir': '/home/mark.miesch/data/lasco_monthly/c3/2014_01/',
    'file': '/17/33385593.fts',
    'bgfile': 'background.fts'
    },
    3: { # CME model 0 from E. Provornikova & A. Malanushenko
    'instrument': 'ModelHAO0',
    'detector': 'original', # the detector field is used to specify the FOV
    'dir': '/home/mark.miesch/data/anny/CME0/pos-30/dcmer_030W_bang_0000_fits/tB/',
    'file': 'frame_0050.fits',
    'bgfile': 'frame_0000.fits'
    }
}

class snapshot:
    """A snapshot is defined as a single coronagraph images for a particular instrument, detector, and time"""

    def __init__(self, file = None, bgfile = None, instrument = None, \
                       detector = None, rescale = False):

        if (file is None):
            print("Incomplete argument list: Using default test case")
            s = testsnap[1]
            instrument = s['instrument']
            detector = s['detector']
            file = s['dir'] + s['file']

        if instrument is None:
            self.instrument = 'None'
        else:
            try:
                self.instrument = instrument.value
            except:
                self.instrument = instrument

        if detector is None:
            self.detector = 'None'
        else:
            try:
                self.detector = detector.value
            except:
                self.detector = detector

        self.file = file
        self.bgfile = bgfile

        # Now read in the data from files.  The assumption here is that there is one image per file.
        # If that is not the case, then we can generalize this as needed.
        try:
            hdu = fits.open(file)[0]
        except FileNotFoundError as fnf_err:
            logging.exception(red+"Fatal Error in CIMP.Event.event constructor: {}".format(fnf_err)+cend)
            raise
        except ValueError as val_err:
            logging.exception(red+"Fatal Error in CIMP.Event.event constructor: Incorrect read API for {}: {}".format(file, val_err)+cend)
            raise
        except BaseException as err:
            logging.exception(red+"Fatal Error in CIMP.Event.event constructor: reading file {} : {}".format(file, err)+cend)
            raise

        self.rawdata = hdu.data.astype('float')
        self.data = self.rawdata.copy()
        if rescale:
            self.rescale()
        self.header = hdu.header

        # adjustments for simulation data
        if 'modelhao' in self.instrument.lower():
            self.header['cunit1'] = 'arcsec'
            self.header['cunit2'] = 'arcsec'

        self.nx = self.header['NAXIS1']
        self.ny = self.header['NAXIS2']

        self.dmap = sunpy.map.Map(self.data, self.header)

        self.time = self.dmap.date

        if self.bgfile is None:
            self.background = None
        else:
            bghdu = fits.open(bgfile)[0]
            self.background = bghdu.data

    @classmethod
    def testcase(cls, case = 1):
        """This is an alternative constructor that creates an Event object based on predefined reference cases.  It takes a single integer as an input, which is the number of the desired test case"""

        s = testsnap[case]
        file = s['dir'] + s['file']

        if s['bgfile'] is None:
            bgfile = None
        else:
            bgfile = s['dir'] + s['bgfile']

        return cls(file, bgfile, s['instrument'], s['detector'])

    def min(self):
        nonzero_pix = np.ma.masked_equal(self.data, 0.0, copy = False)
        return nonzero_pix.min()

    def max(self):
        nonzero_pix = np.ma.masked_equal(self.data, 0.0, copy = False)
        return nonzero_pix.max()

    def min_positive(self):
        positive_pix = np.masked_less_equal(self.data, 0.0, copy = False)
        return positive_pix.min()

    def clip(self, limits, rescale = False):
        self.data = self.data.clip(min = limits[0], max = limits[1])
        if rescale:
            self.rescale()

    def map(self):
        return sunpy.map.Map(self.data, self.header)

    def subtract_background(self, rescale = False, cliprange = 'image'):
        self.data = self.data - self.background

        if rescale:
            self.rescale(cliprange = cliprange)

    def background_ratio(self, rescale = False, cliprange = 'image'):
        """
        This method will compute a ratio relative to the background image and optionally re-scale.  This is similar to NRL's pipeline for creating the daily "pretties" for LASCO, STERO, and PSP/WISPR.
        """

        self.data = np.where(self.background <= 0.0, 0.0, \
                             self.data / self.background)

        if rescale:
            self.rescale(cliprange = cliprange)

    def enhance(self, clip = None, point = 'omr', detail = 'mgn', \
                noise_filter = 'omr', equalize = True):
        """
        Combination of multiple processing steps for plotting
        clip (optional): 2-element tuple specifying the range to clip the data
        point: specify the point filter to use.  Current options are 'omr' or 'median'
        detail: algorithm to use for detail enhancement
        noise_filter: algorithm to use to remove noise after the detail filter has been applied
        """

        print(f"Point Filter: {point}")
        print(f"Detail Enhancement: {detail}")
        print(f"Noise filter: {noise_filter}")

        if point == 'omr':
            a = Enhance.omr(self.data, rescaleim = False)
        elif point == 'median':
            a = Enhance.bright_point_filter(self.data, rescaleim = False)
        else:
            a = self.data

        # contrast stretching via clipping
        if clip is not None:
            a = Enhance.clip(a, min = clip[0], max = clip[1])

        # various techniques to bring out detail
        a = Enhance.detail(a, self.header, filter = detail, \
                           params = [0.8, 1.5])

        # optionally remove noise
        a = Enhance.denoise(a, noise_filter)

        # adaptive equalization
        if detail == 'mgn':
            equalize = False

        if equalize:
            a = Enhance.equalize(a)

        self.data = a

    def equalize(self):
        self.data = Enhance.equalize(self.data)

    def rescale(self, cliprange = 'image'):
        self.data = exposure.rescale_intensity(self.data,
                                               in_range = cliprange, \
                                               out_range = (0,1))

    def point_filter(self, threshold = 2.0, radius = 20):
        self.data = Enhance.bright_point_filter(self.data, \
                    threshold = threshold, radius = radius)
        self.rescale()

    def nrgf(self):
        """
        Apply a normalized radial gradient filter
        """

        amap = Enhance.nrgf(self.map(), self.instrument,
                                        self.detector)
        self.data = amap.data

    def fnrgf(self, order = 20, rmix = [1,15]):
        """
        Apply a Fourier normalized radial gradient filter
        """

        amap = Enhance.fnrgf(self.map(), self.instrument,
                             self.detector, order, rmix)
        self.data = amap.data

    def powerlaw(self, rescale = False):
        self.data = Enhance.powerlaw(self.data)
        if rescale:
            self.rescale()

    def mask_annulus(self, rmin = 0.0, rmax = None):
        Enhance.mask_annulus(self.data, rmin = rmin, rmax = rmax)

    def mask_background(self, rmin = 0.0, rmax = None, nonzero = True):
        # If computing the ratio with the background, you don't want to
        # set the missing value to zero.  But, this does make sense if
        # you're computing a difference

        if nonzero:
            pos_pix = np.ma.masked_less_equal(self.background, 0.0)
            missingval = pos_pix.min()
            self.background = np.where(self.background > 0.0, self.background, missingval)
        else:
            missingval = 0.0

        Enhance.mask_annulus(self.background, rmin = rmin, rmax = rmax,
                             missingval = missingval)

    def downsample(self, rescale = False):
        self.data = Enhance.downsample(self.data)
        self.nx = self.data.shape[0]
        self.ny = self.data.shape[1]
        self.header['NAXIS1'] = self.nx
        self.header['NAXIS2'] = self.ny
        if rescale:
            self.rescale()

    def valid(self, ref = None, tolerance = None, diff_ratio = 100.0):
        """
        This identifies and flags corrupted images by comparing the data to a reference image that is passed as ref (a 2D numpy array).  In a movie sequence, ref would typically be the previous frame.
        Images are flagged as corrupted if more than `diff_ratio` percent of their pixels are different from the reference image, within a specified `tolerance`.  The default diff_ratio is set to be 100%,
        which means that, by default, no frames are declared invalid.  
        The default tolerance is currently set to be 0.5 times the median non-zero value of the reference image.
        """

        if ref is None:
            return True

        if tolerance is None:
            valid_pix = np.ma.masked_equal(ref, 0.0, copy = False)
            tolerance = 0.5*np.nanmedian(np.absolute(valid_pix))

        d = astropy.io.fits.ImageDataDiff(self.data, ref, rtol = tolerance)

        if (100*d.diff_ratio > diff_ratio):
            print(yellow+f"Corrupted image {self.file} {tolerance} {d.diff_ratio*100} {np.min(ref)} {np.max(ref)}" + cend)
            return False
        else:
            return True

    def __str__(self):
        return (f'Instrument = {self.instrument} \n'
               f'Detector   = {self.detector}\n'
               f'Time = {self.time}\n')

    def __repr__(self):
        return (f'Instrument = {self.instrument} \n'
               f'Detector = {self.detector}\n'
               f'file = {self.file}\n'
               f'time = {self.time}')