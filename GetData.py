"""
This is a simple script to explore the contents of a file, including metadata and the header format
"""

import numpy as np
from sunpy.net import attrs as a
from sunpy.net import Fido
import sunpy.io

#------------------------------------------------------------------------------
# define instrument and time range

instrument = a.Instrument.lasco
detector = a.Detector.c2
timerange = a.Time('2013/05/17 9:00:00', '2013/05/17 11:30:00')

#------------------------------------------------------------------------------

dbpath = dir + '/' + instrument.value + '_' + detector.value + '/'

qr = Fido.search(timerange, instrument, detector)
files = Fido.fetch(qr, path = dbpath)

file = files[0]
data, header = sunpy.io.fits.read(file)[0]

