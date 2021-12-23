"""
CIMP Event module  
"""

import astropy.units as u
import datetime
import logging
import numpy as np
import os
import sunkit_image.radial as radial
import sunpy.map
import sunpy.io

from io import RawIOBase
from skimage import exposure
from sunkit_image.utils import equally_spaced_bins
from sunpy.net import Fido
from sunpy.net import attrs as a

# for warning / error statements; print red, yellow text to terminal
red = '\033[91m'
yellow = '\033[93m'
cend = '\033[0m'

testevent = {
    1: {
    'instrument': a.Instrument.lasco,
    'detector': a.Detector.c2,
    'dir': '/home/mark.miesch/sunpy/data/LASCO/',
    'files': ['22605555.fts','22605556.fts','22605557.fts',
              '22605558.fts','22605562.fts','22605563.fts','22605564.fts','22605565.fts',
              '22605566.fts','22605567.fts', '22605568.fts','22605569.fts','22605570.fts']
    },
    2: {
    'instrument': a.Instrument.lasco,
    'detector': a.Detector.c3,
    'dir': '/home/mark.miesch/sunpy/data/LASCO/',
    'files': ['35473922.fts', '35473923.fts']
    }
}

fov = {
    'LASCO-C2' : (1.5, 6.0),
    'LASCO-C3' : (3.7, 30.0)
}

class event:
    """An event is defined as a series of coronagraph images for a particular instrument, detector, and time interval"""

    def __init__(self, instrument = None, detector = None, files = None):

        if (instrument is None or detector is None or files is None):
            print("Incomplete argument list: Using default test case")
            e = testevent[1]
            instrument = e['instrument']
            detector = e['detector']
            files = list(e['dir']+file for file in e['files'])

        self.instrument = instrument
        self.detector = detector
        self._files = files

        # Now read in the data from files.  The assumption here is that there is one image per file.
        # If that is not the case, then we can generalize this as needed.
        self._frames = []
        self.header = []
        self.times = []
        for file in files:
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
               t = header['DATE-OBS'].replace('/','-') + ' ' + header['TIME-OBS']
               time = (datetime.datetime.fromisoformat(t))
            except KeyError as key_err:
                logging.exception(red+"Fatal Error in CIMP.Event.event constructor: header key error {}".format(key_err)+cend)
            except ValueError as val_err:
                logging.exception(red+"Fatal error in CIMP.Event.event constructor: time conversion {}".format(val_err))
            except Exception as e:
                logging.exception(red+'Fatal error in CIMP.Event.event constructor: reading time {}'.format(e)+cend)

            if len(self._frames) == 0:
                self._frames.append(data.astype(float))
            else:
                self._frames.append(data.astype(float) - self._frames[0])

            self.header.append(header)
            self.times.append(time)

        self.nframes = len(self._frames)

    @classmethod
    def testcase(cls, case = 1):
        """This is an alternative constructor that creates an Event object based on predefined reference cases.  It takes a single integer as an input, which is the number of the desired test case"""

        e = testevent[case]
        files = list(e['dir']+file for file in e['files'])

        return cls(e['instrument'], e['detector'], files)

    @classmethod
    def fromtime(cls, instrument = a.Instrument.lasco, detector = a.Detector.c2, 
                 timerange = a.Time('2016/09/06 9:00:00', '2016/09/06 12:00:00'),
                 dir = os.path.expanduser('~')+'/sunpy/data'):
        """This is an alternative constructor that creates an event object based on a selected time interval, specified as a sunpy time range."""

        dbpath = dir + '/' + instrument.value

        qr = Fido.search(timerange, instrument, detector)

        files = Fido.fetch(qr, path = dbpath)

        print(f"CIMP Event constructor downloaded the following files:")
        for file in files:
            print(file)

        return cls(instrument, detector, files)

    def duration(self):
        return (self.times[self.nframes-1] - self.times[0])

    def map(self, idx):
        return sunpy.map.Map(self._frames[idx], self.header[idx])

    def sum(self):
        """
        Long exposure sum of all frames greater than 1.
        Returned as a sunpy map, ready for plotting
        """
        s = self._frames[1]
        for i in np.arange(2,self.nframes):
            s += self._frames[i]
        return sunpy.map.Map(s, self.header[0])

    def enhance(self, clip = None):
        """
        Enhance image frames for plotting
        clip (optional): 2-element tuple specifying the range to clip the data
        """

        # contrast stretching via clipping
        for i in np.arange(1,self.nframes):

            if clip is None:
                a = self._frames[i]
            else:
                a = self._frames[i].clip(min = clip[0], max = clip[1])

            vmin = np.amin(a)
            vmax = np.amax(a)
            im = (a - vmin)/(vmax - vmin)

            # adaptive equalization
            imeq = exposure.equalize_adapthist(im)
            
            self._frames[i] = (vmax - vmin)*imeq.astype(float) + vmin

    def nrgf(self):
        """
        Apply a normalized radial gradient filter to all frames > 0
        """
        
        myfov = fov[self.instrument.value+'-'+self.detector.value]
        
        edges = equally_spaced_bins(myfov[0], myfov[1])
        edges *= u.R_sun

        for i in np.arange(1, self.nframes):
            map = radial.nrgf(self.map(i), edges)
            self._frames[i] = map.data

    def __len__(self):
        return self.nframes

    def __getitem__(self,idx):
        return self._frames[idx]

    def __str__(self):
        return (f'Instrument = {self.instrument.value} \n'
               f'Detector   = {self.detector.value}\n'
               f'Start time = {self.times[0]}\n'
               f'End time   = {self.times[self.nframes-1]}\n'
               f'Duration   = {self.duration()}\n'
               f'Nframes    = {self.nframes}')

    def __repr__(self):
        return (f'Instrument = {self.instrument} \n'
               f'Detector = {self.detector}\n'
               f'Nframes = {self.nframes}\n'
               f'files = {self._files}\n'
               f'times = {self.times}')