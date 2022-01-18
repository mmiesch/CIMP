"""
This script is similar to its counterpart, NoiseRemoval.py, but it operates on individual images as opposed to difference images.
"""

import copy
import numpy as np
from CIMP import Snapshot as snap
import sunpy.map
from sunpy.net import attrs as a
import matplotlib.pyplot as plt
import astroscrappy
import noisegate as ng

from skimage import exposure
from skimage.filters import median
from skimage.filters.rank import enhance_contrast
from skimage.morphology import disk, remove_small_objects, white_tophat
from skimage.restoration import (denoise_tv_chambolle, denoise_bilateral,
                                 denoise_wavelet, estimate_sigma, denoise_nl_means)

import scipy.ndimage

def remove_outliers(im, radius = 2, threshold = 50):
    medim = median(im,disk(radius))
    outliers = ( (im > medim + threshold) |
                 (im < medim - threshold) )
    out = np.where(outliers, medim, im)
    return out


plotcase = 1

if plotcase == 1:
    instrument = a.Instrument.lasco
    detector = a.Detector.c3
    file = '/home/mark.miesch/sunpy/data/lasco_c3/32305543.fts'
    bgfile = '/home/mark.miesch/data/lasco_ssw/3m_clcl_120716.fts'
    nrgf = False
    scale = None

else:
    print("specify a valid plotcase")
    exit()    


#======================================================================
# get image and remove background

x = snap.snapshot(instrument, detector, file)
a = copy.deepcopy(x.data.astype('float'))
asc = exposure.rescale_intensity(a)

x.background_normalize(bgfile)
b = copy.deepcopy(x.data.astype('float'))

#if nrgf:
#    x.nrgf()

# for experimenting
#timerange = a.Time('2016/09/06 8:00:00', '2016/09/06 12:00:00')
#x = ev.event.fromtime(a.Instrument.lasco, a.Detector.c2, timerange)

print(80*'-')
print(x)
print(80*'-')

#======================================================================
#======================================================================
# plot

fig = plt.figure(figsize=[24,12])

p = exposure.equalize_adapthist(asc)
pmap = sunpy.map.Map(p, x.header)
ax = fig.add_subplot(2,3,1,projection=pmap)
pmap.plot(vmin = 0.0, vmax = 1.0)

psc = exposure.rescale_intensity(x.background.astype('float'))
p = exposure.equalize_adapthist(psc)
pmap = sunpy.map.Map(p, x.header)
ax = fig.add_subplot(2,3,2,projection=pmap)
pmap.plot(vmin = 0.0, vmax = 1.0)

p = exposure.equalize_adapthist(b)
pmap = sunpy.map.Map(p, x.header)
ax = fig.add_subplot(2,3,3,projection=pmap)
pmap.plot(vmin = 0.0, vmax = 0.1)

#======================================================================
plt.show()