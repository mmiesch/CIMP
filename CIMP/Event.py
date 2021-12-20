"""
CIMP Event module  
"""

from sunpy.net import Fido
from sunpy.net import attrs as a

testevent = {
    1: {
    'instrument': a.Instrument.lasco,
    'detector': a.Detector.c2,
    'dir': '/home/mark.miesch/sunpy/data/lasco/',
    'files': ['22605567.fts', '22605568.fts']
    },
    2: {
    'instrument': a.Instrument.lasco,
    'detector': a.Detector.c3,
    'dir': '/home/mark.miesch/sunpy/data/lasco/',
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

        self.nframes = 0
        self.frames = None

    @classmethod
    def testcase(cls, case = 1):
        """This is an alternative constructor that creates an Event object based on predefined reference cases.  It takes a single integer as an input, which is the number of the desired test case"""

        e = testevent[case]
        files = list(e['dir']+file for file in e['files'])

        return cls(e['instrument'], e['detector'], files)

    @classmethod
    def fromtime(cls, instrument = a.Instrument.lasco, detector = a.Detector.c2, 
                 timerange = a.Time('2016/09/06 10:00:00', '2016/09/06 11:00:00')):
        """This is an alternative constructor that creates an event object based on a selected time interval, specified as a sunpy time range."""

        # placeholder for testing
        instrument = a.Instrument.lasco
        detector = a.Detector.c2
        dir = '/home/mark.miesch/sunpy/data/lasco/'
        filelist = ['22605567.fts', '22605568.fts']
        files = list(dir+file for file in filelist)

        return cls(instrument, detector, files)

    def __str__(self):
        return (f'Instrument = {self.instrument.value} \n'
               f'Detector = {self.detector.value}\n')

    def __repr__(self):
        return (f'Instrument = {self.instrument} \n'
               f'Detector = {self.detector}\n'
               f'files = {self._files}')