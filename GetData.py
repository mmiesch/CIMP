"""
This is a simple script to explore the contents of a file, including metadata and the header format
"""

import numpy as np
from sunpy.net import attrs as a
from sunpy.net import Fido
import sunpy.io
import sunpy.map

#------------------------------------------------------------------------------
# define instrument and time range

dcase = 6

if dcase == 1:

   # usually it's not necessary to specify the source.  
   # An exception is STEREO, where you have to specify the A or B spacecreaft
   source = None
   instrument = a.Instrument.lasco
   detector = a.Detector.c2
   timerange = a.Time('2013/05/17 9:00:00', '2013/05/17 11:00:00')

elif dcase == 2:
   source = None
   instrument = a.Instrument.lasco
   detector = a.Detector.c3
   timerange = a.Time('2013/05/17 9:30:00', '2013/05/17 11:00:00')

elif dcase == 3:
   source = "STEREO_A"
   instrument = a.Instrument.secchi
   detector = a.Detector.cor1
   timerange = a.Time('2013/05/17 9:00:00', '2013/05/17 13:30:00')

elif dcase == 4:
   source = "STEREO_B"
   instrument = a.Instrument.secchi
   detector = a.Detector.cor1
   timerange = a.Time('2013/05/17 9:00:00', '2013/05/17 13:30:00')

elif dcase == 5:
   source = "STEREO_A"
   instrument = a.Instrument.secchi
   detector = a.Detector.cor2
   timerange = a.Time('2013/05/17 9:00:00', '2013/05/17 13:30:00')

elif dcase == 6:
   source = "STEREO_B"
   instrument = a.Instrument.secchi
   detector = a.Detector.cor2
   timerange = a.Time('2013/05/17 9:00:00', '2013/05/17 13:30:00')

else:
    print("specify a valid dcase")
    exit(1)


#------------------------------------------------------------------------------
basedir = '/home/mark.miesch/sunpy/data/'

dbpath = basedir + instrument.value.lower() + '_' + detector.value.lower() + '/'

if source is None:
    qr = Fido.search(timerange, instrument, detector)
else:
    qr = Fido.search(timerange, a.Source(source), instrument, detector)


files = Fido.fetch(qr, path = dbpath)
files.sort()

print(files)

file = files[0]
data, header = sunpy.io.fits.read(file)[0]

dmap = sunpy.map.Map(data, header)

if source is None:
    print(f"{instrument.value} / {detector.value} : {dmap.date.value}")
else:
    print(f"{source} / {instrument.value} / {detector.value} : {dmap.date.value}")
