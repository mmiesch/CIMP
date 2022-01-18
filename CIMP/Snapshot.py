"""
CIMP Event module  
"""

import astropy.units as u
import datetime
import logging
from astropy.utils import data_info
import numpy as np
import os
import sunkit_image.radial as radial
import sunpy.map
import sunpy.io

from io import RawIOBase
from skimage import exposure
from skimage.filters.rank import median, enhance_contrast_percentile
from skimage.morphology import disk
from skimage.restoration import (denoise_tv_chambolle, denoise_nl_means)
from sunkit_image.utils import equally_spaced_bins
from sunpy.net import Fido
from sunpy.net import attrs as a

# for warning / error statements; print red, yellow text to terminal
red = '\033[91m'
yellow = '\033[93m'
cend = '\033[0m'

testsnap = {
    1: {
    'instrument': a.Instrument.lasco,
    'detector': a.Detector.c2,
    'dir': '/home/mark.miesch/sunpy/data/LASCO/',
    'file': '22605555.fts'
    },
    2: {
    'instrument': a.Instrument.lasco,
    'detector': a.Detector.c3,
    'dir': '/home/mark.miesch/sunpy/data/LASCO/',
    'file': '32473914.fts'
    }
}

fov = {
    'LASCO-C2' : (1.5, 6.0),
    'LASCO-C3' : (3.7, 30.0)
}

class snapshot:
    """A snapshot is defined as a single coronagraph images for a particular instrument, detector, and time"""

    def __init__(self, instrument = None, detector = None, file = None):

        if (instrument is None or detector is None or file is None):
            print("Incomplete argument list: Using default test case")
            s = testsnap[1]
            instrument = s['instrument']
            detector = s['detector']
            file = s['dir'] + s['file']

        self.instrument = instrument
        self.detector = detector
        self.file = file

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

        try:
            print(f"HEADER TIME {file} {header['DATE-OBS']}")
            #if (self.instrument.value == 'LASCO' and self.detector.value == 'C2'): 
            t = header['DATE-OBS'].replace('/','-') + ' ' + header['TIME-OBS']
            time = datetime.datetime.fromisoformat(t)
            #else:
            #    """
            #    Careful here - it looks like the time stamp on C3 files can be different - like C2 or not
            #    """
            #    t = header['DATE-OBS'].replace('T',' ')
            #    time = datetime.datetime.fromisoformat(t)
        except KeyError as key_err:
            logging.exception(red+"Fatal Error in CIMP.Event.event constructor: header key error {}".format(key_err)+cend)
        except ValueError as val_err:
            logging.exception(red+"Fatal error in CIMP.Event.event constructor: time conversion {} : {}".format(t, val_err)+cend)
        except Exception as e:
            logging.exception(red+'Fatal error in CIMP.Event.event constructor: reading time {}'.format(e)+cend)

        self.data = data.astype('float')
        self.header = header 
        self.time = time

    @classmethod
    def testcase(cls, case = 1):
        """This is an alternative constructor that creates an Event object based on predefined reference cases.  It takes a single integer as an input, which is the number of the desired test case"""

        s = testsnap[case]
        file = s['dir'] + s['file']

        return cls(s['instrument'], s['detector'], file)

    def map(self):
        return sunpy.map.Map(self.data, self.header)

    def background_normalize(self, bgfile = None):
        """
        Given a filename that contains a background image, this method will compute a ratio relative to the background image and re-scale.  This is similar to NRL's pipeline for creating the daily "pretties" for LASCO, STERO, and PSP/WISPR.
        """
        
        # Currently a filename is a required parameter
        if (bgfile is None):
            print("Error in Snapshot.background_normalize(): you must specify a filename that contains the background image")
            exit()
        
        self.background, bheader = sunpy.io.fits.read(bgfile)[0]

        rat = np.where(self.background <= 0.0, 0.0, self.data/self.background)
        self.data = exposure.rescale_intensity(rat)


    def enhance(self, clip = None, noise_removal = 'tv'):
        """
        Enhance image frames for plotting
        clip (optional): 2-element tuple specifying the range to clip the data
        """
        vmin = clip[0]
        vmax = clip[1]

        # contrast stretching via clipping

        if clip is None:
            a = self.data
        else:
            a = self.data.clip(min = vmin, max = vmax)

        im = (a - vmin)/(vmax - vmin)

        # remove noise
        if noise_removal == 'tv':
            imdn = denoise_tv_chambolle(im, weight = 0.2)
        elif noise_removal == 'mediam':
            imdn = median(disk(1))
        else:
            imdn = denoise_nl_means(im, patch_size = 4)

        imdn = (imdn - np.amin(imdn))/(np.amax(imdn) - np.amin(imdn))
        imc = enhance_contrast_percentile(imdn, disk(2), p0=.1, p1=.9)

        # adaptive equalization
        imeq = exposure.equalize_adapthist(imc)
        
        self.data = (vmax - vmin)*imeq + vmin

    def nrgf(self):
        """
        Apply a normalized radial gradient filter to all frames > 0
        """
        
        myfov = fov[self.instrument.value+'-'+self.detector.value]
        
        edges = equally_spaced_bins(myfov[0], myfov[1])
        edges *= u.R_sun

        map = radial.nrgf(self.map(), edges)
        self.data = map.data

    def __str__(self):
        return (f'Instrument = {self.instrument.value} \n'
               f'Detector   = {self.detector.value}\n'
               f'Time = {self.time}\n')

    def __repr__(self):
        return (f'Instrument = {self.instrument} \n'
               f'Detector = {self.detector}\n'
               f'file = {self.file}\n'
               f'time = {self.time}')