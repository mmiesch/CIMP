import numpy as np
from sunpy.net import attrs as a
from sunpy.net import Fido
import sunpy.io
import sunpy.map

dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04'

instrument = a.Instrument.lasco
detector = a.Detector.c3
level = a.Level.two

month = '2012/04'
Ndays = 30

for d in np.arange(1,Ndays+1):
    day = f'{d:02d}'
    date = month+'/'+day
    t1 = date + ' 0:00:00'
    t2 = date + ' 11:59:59'
    timerange = a.Time(t1, t2)
    qr = Fido.search(timerange, instrument, detector, level)
    Fido.fetch(qr,path=dir+'/'+day)

