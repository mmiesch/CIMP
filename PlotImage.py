
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import sunpy.io

#dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/'
#file = dir+'03/32295364.fts'

dir = '/home/mark.miesch/data/anny/CME0/pos-30/dcmer_030W_bang_0000_fits/tB/'
file = dir+'frame_0050.fits'

data, header = sunpy.io.fits.read(file)[0]

try:
    dmap = sunpy.map.Map(data, header)
    dmap.plot()
except:
    plt.imshow(data)

plt.show()