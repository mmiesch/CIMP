"""
CIMP Event module
"""

import astropy.units as u
import datetime
import logging
from astropy.utils import data_info
import numpy as np
import os
import sunpy.map
import sunpy.io

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
    }
}

class snapshot:
    """A snapshot is defined as a single coronagraph images for a particular instrument, detector, and time"""

    def __init__(self, file = None, bgfile = None, instrument = None, \
                       detector = None):

        if (file is None):
            print("Incomplete argument list: Using default test case")
            s = testsnap[1]
            instrument = s['instrument']
            detector = s['detector']
            file = s['dir'] + s['file']

        if instrument is None:
            self.instrument = 'None'
        else:
            self.instrument = instrument.value

        if detector is None:
            self.detector = 'None'
        else:
            self.detector = detector.value

        self.file = file
        self.bgfile = bgfile

        # Now read in the data from files.  The assumption here is that there is one image per file.
        # If that is not the case, then we can generalize this as needed.
        try:
            data, header = sunpy.io.fits.read(file)[0]
        except FileNotFoundError as fnf_err:
            logging.exception(red+"Fatal Error in CIMP.Event.event constructor: {}".format(fnf_err)+cend)
            raise
        except ValueError as val_err:
            logging.exception(red+"Fatal Error in CIMP.Event.event constructor: Incorrect read API for {}: {}".format(file, val_err)+cend)
            raise
        except BaseException as err:
            logging.exception(red+"Fatal Error in CIMP.Event.event constructor: reading file {} : {}".format(file, err)+cend)
            raise

        self.rawdata = data.astype('float')
        self.data = exposure.rescale_intensity(self.rawdata, out_range=(0,1))
        self.header = header

        # adjustments for simulation data
        if not 'cunit1' in self.header.keys():
            self.header['cunit1'] = 'arcsec'
            self.header['cunit2'] = 'arcsec'

        self.dmap = sunpy.map.Map(self.data, self.header)

        self.time = self.dmap.date

        if self.bgfile is None:
            self.background = None
        else:
            self.background, bgheader = sunpy.io.fits.read(bgfile)[0]

    @classmethod
    def testcase(cls, case = 1):
        """This is an alternative constructor that creates an Event object based on predefined reference cases.  It takes a single integer as an input, which is the number of the desired test case"""

        s = testsnap[case]
        file = s['dir'] + s['file']
        bgfile = s['dir'] + s['bgfile']

        return cls(file, bgfile, s['instrument'], s['detector'])

    def clip(self,limits):
        self.data = self.data.clip(min = limits[0], max = limits[1])
        self.rescale()

    def map(self):
        return sunpy.map.Map(self.data, self.header)

    def subtract_background(self):
        self.data = self.rawdata - self.background
        self.rescale()

    def background_normalize(self):
        """
        This method will compute a ratio relative to the background image and re-scale.  This is similar to NRL's pipeline for creating the daily "pretties" for LASCO, STERO, and PSP/WISPR.
        """

        self.data = np.where(self.background <= 0.0, 0.0, \
                             self.rawdata / self.background)

        self.rescale()

    def enhance(self, clip = None, detail = 'mgn', noise_filter = 'bregman'):
        """
        Combination of multiple processing steps for plotting
        clip (optional): 2-element tuple specifying the range to clip the data
        """

        print(f"Detail Enhancement: {detail}")
        print(f"Noise filter: {noise_filter}")

        a = Enhance.bright_point_filter(self.data, rescaleim = False)

        # contrast stretching via clipping
        if clip is not None:
            a = Enhance.clip(a, min = clip[0], max = clip[1])

        # various techniques to bring out detail
        a = Enhance.detail(a, self.header, filter = 'mgn', \
                           params = [0.8, 1.5])

        # optionally remove noise
        a = Enhance.denoise(a, noise_filter)

        # adaptive equalization
        if (detail != 'mgn'):
            a = Enhance.equalize(a)

        self.data = a

    def equalize(self):
        self.data = Enhance.equalize(self.data)

    def rescale(self):
        self.data = exposure.rescale_intensity(self.data, out_range=(0,1))

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

    def mask_annulus(self, rmin = 0.0, rmax = None):
        Enhance.mask_annulus(self.data, rmin = rmin, rmax = rmax)

    def __str__(self):
        return (f'Instrument = {self.instrument} \n'
               f'Detector   = {self.detector}\n'
               f'Time = {self.time}\n')

    def __repr__(self):
        return (f'Instrument = {self.instrument} \n'
               f'Detector = {self.detector}\n'
               f'file = {self.file}\n'
               f'time = {self.time}')