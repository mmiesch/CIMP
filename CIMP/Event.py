"""
CIMP Event module  
"""

from io import RawIOBase
import logging
import os
from sunpy.net import Fido
from sunpy.net import attrs as a
import sunpy.io

# for warning / error statements; print red, yellow text to terminal
red = '\033[91m'
yellow = '\033[93m'
cend = '\033[0m'

testevent = {
    1: {
    'instrument': a.Instrument.lasco,
    'detector': a.Detector.c2,
    'dir': '/home/mark.miesch/sunpy/data/LASCO/',
    'files': ['22605566.fts','22605567.fts', '22605568.fts','22605569.fts','22605570.fts']
    },
    2: {
    'instrument': a.Instrument.lasco,
    'detector': a.Detector.c3,
    'dir': '/home/mark.miesch/sunpy/data/LASCO/',
    'files': ['35473922.fts', '35473923.fts']
    }
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
        self.nframes = 0
        self.frames = []
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
               t = header['DATE'].replace('/','-')
               time = (datetime.datetime.fromisoformat(t))
            except Exception as e:
                logging.exception('Unexpected error in CIMP.Event.event constructor: reading time')

            self.frames.append(data.astyp(float64))
            self.times.append(time)
            self.nframes += 1

    @classmethod
    def testcase(cls, case = 1):
        """This is an alternative constructor that creates an Event object based on predefined reference cases.  It takes a single integer as an input, which is the number of the desired test case"""

        e = testevent[case]
        files = list(e['dir']+file for file in e['files'])

        return cls(e['instrument'], e['detector'], files)

    @classmethod
    def fromtime(cls, instrument = a.Instrument.lasco, detector = a.Detector.c2, 
                 timerange = a.Time('2016/09/06 10:00:00', '2016/09/06 11:00:00'),
                 dir = os.path.expanduser('~')+'/sunpy/data'):
        """This is an alternative constructor that creates an event object based on a selected time interval, specified as a sunpy time range."""

        dbpath = dir + '/' + instrument.value

        qr = Fido.search(timerange, instrument, detector)

        files = Fido.fetch(qr, path = dbpath)

        print(f"CIMP Event constructor downloaded the following files:")
        for file in files:
            print(file)

        return cls(instrument, detector, files)

    def __str__(self):
        return (f'Instrument = {self.instrument.value} \n'
               f'Detector = {self.detector.value}\n'
               f'Start time = {self.times[0]}')

    def __repr__(self):
        return (f'Instrument = {self.instrument} \n'
               f'Detector = {self.detector}\n'
               f'files = {self._files}'
               f'times = {self.times}')