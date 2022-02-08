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
import noisegate as ng
import sunkit_image.enhance

from io import RawIOBase
from skimage import exposure
from skimage.filters import median
from skimage.filters.rank import enhance_contrast_percentile
from skimage.morphology import disk
from skimage.restoration import (denoise_tv_chambolle, denoise_tv_bregman, 
                                 denoise_nl_means)
from sunkit_image.utils import equally_spaced_bins
from sunpy.net import Fido
from sunpy.net import attrs as a

# for warning / error statements; print red, yellow text to terminal
red = '\033[91m'
yellow = '\033[93m'
cend = '\033[0m'

testevent = {
    1: {
    'source': None,
    'instrument': a.Instrument.lasco,
    'detector': a.Detector.c2,
    'dir': '/home/mark.miesch/sunpy/data/lasco_c2/',
    'files': list('226055'+num+'.fts' for num in ['55','56','57','58',
      '62','63','64','65','66','67','68','69','70','71'])
    },
    2: {
    'source': None,
    'instrument': a.Instrument.lasco,
    'detector': a.Detector.c3,
    'dir': '/home/mark.miesch/sunpy/data/lasco_c3/',
    'files': list(str(num)+'.fts' for num in np.arange(32473914, 32473944))
    },
    3: {
    'source': None,
    'instrument': a.Instrument.lasco,
    'detector': a.Detector.c2,
    'dir': '/home/mark.miesch/sunpy/data/lasco_c2/',
    'files': list(str(num)+'.fts' for num in np.arange(22459503, 22459510))
    },
    4: {
    'source': None,
    'instrument': a.Instrument.lasco,
    'detector': a.Detector.c3,
    'dir': '/home/mark.miesch/sunpy/data/lasco_c3/',
    'files': list(str(num)+'.fts' for num in np.arange(32339573, 32339592))
    },
    5: {
    'source': "STEREO_A",
    'instrument': a.Instrument.secchi,
    'detector': a.Detector.cor1,
    'dir': '/home/mark.miesch/sunpy/data/secchi_cor1/',
    'files': list("20130517_"+t+"00_s4c1a.fts" for t in 
             ["0910","0915","0920","0925","0930","0935","0940","0945","0950","0955",
              "1010","1015","1020","1025","1030","1035","1040","1045","1050","1055"])
    },
    6: {
    'source': "STEREO_A",
    'instrument': a.Instrument.secchi,
    'detector': a.Detector.cor2,
    'dir': '/home/mark.miesch/sunpy/data/secchi_cor2/',
    'files': list("20130517_"+t+"00_d4c2a.fts" for t in 
             ["0924","0954","1024","1054","1124","1154"])
    },
    7: {
    'source': "STEREO_B",
    'instrument': a.Instrument.secchi,
    'detector': a.Detector.cor1,
    'dir': '/home/mark.miesch/sunpy/data/secchi_cor1/',
    'files': list("20130517_"+t+"00_s4c1b.fts" for t in 
             ["0910","0915","0920","0925","0930","0935","0940","0945","0950","0955",
              "1010","1015","1020","1025","1030","1035","1040","1045","1050","1055"])
    },
    8: {
    'source': "STEREO_B",
    'instrument': a.Instrument.secchi,
    'detector': a.Detector.cor2,
    'dir': '/home/mark.miesch/sunpy/data/secchi_cor2/',
    'files': list("20130517_"+t+"00_d4c2b.fts" for t in 
             ["0924","0954","1024","1054","1124","1154"])
    }
}

fov = {
    'lasco-c2'   : (1.5,  6.0),
    'lasco-c3'   : (3.7, 30.0),
    'secchi-cor1': (1.5,  4.0),
    'secchi-cor2': (3.0, 15.0)
}

def point_filter(im, threshold = 2.0, radius = 20):
    amag = np.absolute(im)
    amag -= np.min(amag)
    amed = median(amag, disk(radius))
    rob = amag < (1.0+threshold) * amed
    p = im * rob.astype('float')
    return np.where(amag > 0, p, 0.0)

class event:
    """An event is defined as a series of coronagraph images for a particular instrument, detector, and time interval"""

    def __init__(self, instrument = None, detector = None, files = None, source = None):

        if (instrument is None or detector is None or files is None):
            print("Incomplete argument list: Using default test case")
            e = testevent[1]
            instrument = e['instrument']
            detector = e['detector']
            files = list(e['dir']+file for file in e['files'])

        if source is None:
            self.instrument = instrument.value
        else:
            self.instrument = source + '/' + instrument.value
        self.detector = detector.value
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

            if len(self._frames) == 0:
                self._frames.append(data.astype(float))
            else:
                self._frames.append(data.astype(float) - self._frames[0])

            self.header.append(header)

        self.nframes = len(self._frames)

        # different instruments have different headers.  Let Sunpy sort it out with maps
        for i in np.arange(0, self.nframes):
            m = self.map(i)
            try: 
                time = datetime.datetime.fromisoformat(m.date.value)
            except Exception as e:
                logging.exception(red+'Error in CIMP.Event.event constructor : {} : reading time {}'.format(self._file[i], e)+cend)
                time = datetime.datetime(datetime.MAXYEAR,1,1)
            self.times.append(time)

    @classmethod
    def testcase(cls, case = 1):
        """This is an alternative constructor that creates an Event object based on predefined reference cases.  It takes a single integer as an input, which is the number of the desired test case"""

        e = testevent[case]
        files = list(e['dir']+file for file in e['files'])

        return cls(e['instrument'], e['detector'], files, e['source'])

    @classmethod
    def fromtime(cls, instrument = a.Instrument.lasco, detector = a.Detector.c2, 
                 timerange = a.Time('2016/09/06 9:00:00', '2016/09/06 12:00:00'),
                 source = None, dir = os.path.expanduser('~')+'/sunpy/data'):
        """This is an alternative constructor that creates an event object based on a selected time interval, specified as a sunpy time range."""

        dbpath = dir + '/' + instrument.value.lower() + '_' + detector.value.lower() + '/'

        if (instrument.value.lower() == 'secchi' and source is None):
            source = 'STEREO_A'

        if source is None:
            qr = Fido.search(timerange, instrument, detector)
        else:
            qr = Fido.search(timerange, a.Source(source), instrument, detector)


        files = Fido.fetch(qr, path = dbpath)
        files.sort()

        print(f"CIMP Event constructor downloaded the following files:")
        for file in files:
            print(file)

        return cls(instrument, detector, files, source)

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

    def enhance(self, clip = None, noise_filter = 'bregman', detail = 'mgn',
                rmix = [1,15]):
        """
        Enhance image frames for plotting
        clip (optional): 2-element tuple specifying the range to clip the data
        """

        print(f"Detail Enhancement: {detail}")
        print(f"Noise filter: {noise_filter}")

        for i in np.arange(1,self.nframes):

            # filter out bright points
            a = point_filter(self._frames[i])

            # contrast stretching via clipping
            if clip is not None:
                a = a.clip(min = clip[0], max = clip[1])

            # various techniques to bring out detail
            if detail == 'mgn':
                """
                Multiscale Gaussian Noise filter (Morgan & Druckmuller 2014)
                """
                b = sunkit_image.enhance.mgn(a, h = 0.7, gamma = 1.5)
            elif detail == 'fnrgf':
                """
                Fourier Normaling Radial Gradient Filter (Druckmullerova et al 2011)
                """
                myfov = fov[self.instrument.lower()+'-'+self.detector.lower()]
                
                edges = equally_spaced_bins(myfov[0], myfov[1])
                edges *= u.R_sun

                order = 20
                coefs = radial.set_attenuation_coefficients(order)
                amap = sunpy.map.Map(a, self.header[i])
                bmap = radial.fnrgf(amap, edges, order, coefs, ratio_mix = rmix)
                b = bmap.data.clip(min=0.0)
            elif detail == 'contrast':
                asc = exposure.rescale_intensity(a)
                b = enhance_contrast_percentile(asc, disk(2), p0=.1, p1=.9)
            else:
                b = a

            ## optionally remove noise
            if noise_filter == 'tv':
                c = denoise_tv_chambolle(b, weight = 0.2)
            elif noise_filter == 'bregman':
                c = denoise_tv_bregman(b)
            elif noise_removal == 'median':
                c = median(a,disk(1))
            elif noise_filter == 'nl_means"':
                c = denoise_nl_means(b, patch_size = 4)
            else:
                c = b

            # adaptive equalization
            if (detail != 'mgn'):
                csc = exposure.rescale_intensity(c, out_range=(0,1))
                ceq = exposure.equalize_adapthist(csc)
                c = exposure.rescale_intensity(ceq, out_range=(0,1))
            
            self._frames[i] = c

    def nrgf(self):
        """
        Apply a normalized radial gradient filter to all frames > 0
        """
        
        myfov = fov[self.instrument.lower()+'-'+self.detector.lower()]
        
        edges = equally_spaced_bins(myfov[0], myfov[1])
        edges *= u.R_sun

        for i in np.arange(1, self.nframes):
            map = radial.nrgf(self.map(i), edges)
            self._frames[i] = map.data

    def noise_gate(self, cubesize = (3, 12, 12), model = 'hybrid',
                   factor = 2.0):
        """
        Noise Gate filter from DeForest, C.E. 2017, ApJ, 838:155 (10pp)
        """

        print(f"Applying noise gate filter: {cubesize} {model} {factor}")
        dcube = np.zeros((self.nframes, self._frames[0].shape[0], 
                          self._frames[0].shape[1]))
    
        for i in np.arange(1, self.nframes):
            dcube[i-1,:,:] = self._frames[i]

        ng.noise_gate_batch(dcube, cubesize=cubesize, model=model, 
                            factor = factor)

        for i in np.arange(1, self.nframes):
            self._frames[i] = dcube[i-1,:,:]

    def __len__(self):
        return self.nframes

    def __getitem__(self,idx):
        return self._frames[idx]

    def __str__(self):
        return (f'Instrument = {self.instrument} \n'
               f'Detector   = {self.detector}\n'
               f'Resolution = {self._frames[0].shape}\n'
               f'Start time = {self.times[0]}\n'
               f'End time   = {self.times[self.nframes-1]}\n'
               f'Duration   = {self.duration()}\n'
               f'Nframes    = {self.nframes}')

    def __repr__(self):
        return (f'Instrument = {self.instrument} \n'
               f'Detector = {self.detector}\n'
               f'Nframes = {self.nframes}\n'
               f'Resolution = {self._frames[0].shape}\n'
               f'files = {self._files}\n'
               f'times = {self.times}')