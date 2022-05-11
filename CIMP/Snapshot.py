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

class snapshot:
    """A snapshot is defined as a single coronagraph images for a particular instrument, detector, and time"""

    def __init__(self, instrument = None, detector = None, file = None):

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
            logging.exception(red+"Warning in CIMP.Event.event constructor: header key error {}".format(key_err)+cend)
            time = datetime.datetime.now()
        except ValueError as val_err:
            logging.exception(red+"Warning in CIMP.Event.event constructor: time conversion {} : {}".format(t, val_err)+cend)
        except Exception as e:
            logging.exception(red+'Warning in CIMP.Event.event constructor: reading time {}'.format(e)+cend)

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
        if not 'cunit1' in self.header.keys():
            self.header['cunit1'] = 'arcsec'
            self.header['cunit2'] = 'arcsec'

        return sunpy.map.Map(self.data, self.header)

    def background_normalize(self, bgfile = None):
        """
        Given a filename that contains a background image, this method will compute a ratio relative to the background image and re-scale.  This is similar to NRL's pipeline for creating the daily "pretties" for LASCO, STERO, and PSP/WISPR.
        """
        
        # Currently a filename is a required parameter
        if (bgfile is None):
            print("Error in Snapshot.background_normalize(): you must specify a filename that contains the background image")
            exit()
        
        bkg, bheader = sunpy.io.fits.read(bgfile)[0]

        self.background = bkg - np.min(bkg)
        a = self.data - np.min(self.data)
        rat = np.where(self.background <= 0.0, 0.0, a/self.background)
        self.data = exposure.rescale_intensity(rat)

    def enhance(self, clip = None, noise_filter = 'bregman', detail = 'mgn',
                rmix = [1,15]):
        """
        Enhance image frames for plotting
        clip (optional): 2-element tuple specifying the range to clip the data
        """

        print(f"Detail Enhancement: {detail}")
        print(f"Noise filter: {noise_filter}")

        for i in np.arange(1,self.nframes):

            image = Enhance.enhance(self.data,
                                    instrument = self.instrument,
                                    detector = self.detector, 
                                    clip = clip, 
                                    noise_filter = noise_filter, 
                                    detail = detail, 
                                    rmix = rmix)
            self.data = image

    def nrgf(self):
        """
        Apply a normalized radial gradient filter
        """
        
        amap = Enhance.nrgf(self.map(), self.instrument.value, 
                                        self.detector.value)
        self.data = amap.data

    def fnrgf(self, order = 20, rmix = [1,15]):
        """
        Apply a Fourier normalized radial gradient filter
        """
        
        amap = Enhance.fnrgf(self.map(), self.instrument.value, 
                             self.detector.value, order, rmix)
        self.data = amap.data

    def __str__(self):
        return (f'Instrument = {self.instrument} \n'
               f'Detector   = {self.detector}\n'
               f'Time = {self.time}\n')

    def __repr__(self):
        return (f'Instrument = {self.instrument} \n'
               f'Detector = {self.detector}\n'
               f'file = {self.file}\n'
               f'time = {self.time}')