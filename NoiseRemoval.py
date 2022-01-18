"""
The purpose of this script is to start playing around with different python image processing tools
"""

import numpy as np
from CIMP import Event as ev
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
    testcase = 4
    nrgf = False
    idx = 9
    scale = (0.0, 1000.0)

elif plotcase == 2:
    testcase = 4
    nrgf = True
    idx = 3
    scale = (0.0, 4.0)

else:
    print("specify a valid plotcase")
    exit()    


x = ev.event.testcase(testcase)

if nrgf:
    x.nrgf()

# for experimenting
#timerange = a.Time('2016/09/06 8:00:00', '2016/09/06 12:00:00')
#x = ev.event.fromtime(a.Instrument.lasco, a.Detector.c2, timerange)

print(80*'-')
print(x)
print(80*'-')
print(repr(x))
print(80*'-')

# ===================
# pick an image to work with from an event

a = x[idx]

# image as a sunpy map
amap = x.map(idx)

amin = np.amin(a)
amax = np.amax(a)

print(f"Original image range: {amin} to {amax}")

# clipped image
aclip = a.clip(min=scale[0], max=scale[1])

# scaled byte image
#asc = (a - amin)/(amax - amin)
#asc = (255*aclip/np.amax(aclip)).astype('uint8')
asc = aclip/np.amax(aclip)

#======================================================================
# estimate noise

# Estimate the average noise standard deviation across color channels.
sigma_est = estimate_sigma(a, channel_axis=-1, average_sigmas=True)

# Due to clipping in random_noise, the estimate will be a bit smaller than the
# specified sigma.
print(f'Estimated Gaussian noise standard deviation = {sigma_est}')

sigma_est = estimate_sigma(a, channel_axis=-1, average_sigmas=True)
print(f'Estimated Gaussian noise standard deviation (clipped) = {sigma_est}')

#======================================================================
# noise removal

pims = []
titles = []

# total variation filter
#titles.append('tv')
##psc = denoise_tv_chambolle(asc, weight=0.1)
#psc = denoise_tv_chambolle(asc, weight=0.2)
#p = (psc - np.amin(psc))/(np.amax(psc) - np.amin(psc))
#pims.append(p)

# bilateral filter
#titles.append('bilateral')
##psc = denoise_bilateral(asc, sigma_spatial=15)
#psc = denoise_bilateral(asc, sigma_spatial=25)
#p = (psc - np.amin(psc))/(np.amax(psc) - np.amin(psc))
#pims.append(p)

# nlmeans
#titles.append('nlmeans')
#psc = denoise_nl_means(asc, patch_size=4)
#p = (psc - np.amin(psc))/(np.amax(psc) - np.amin(psc))
#pims.append(p)

titles.append('median')
psc = median(asc, disk(1))
p = exposure.rescale_intensity(psc)
pims.append(p)

# wavelet filter
#titles.append('wavelet')
##psc = denoise_wavelet(asc, rescale_sigma=True)
#psc = denoise_wavelet(asc, sigma=0.1, mode='soft',wavelet='haar')
#p = (psc - np.amin(psc))/(np.amax(psc) - np.amin(psc))
#pims.append(p)

#titles.append('medfilt2d')
#psc = scipy.signal.medfilt2d(asc, kernel_size=3)
#p = (psc - np.amin(psc))/(np.amax(psc) - np.amin(psc))
#pims.append(p)#

## Astroscrappy cosmic ray removal
#titles.append('astroscrappy')
#mask, psc = astroscrappy.detect_cosmics(asc, sigclip=2, objlim=2, readnoise=4, verbose=True)
#p = (psc - np.amin(psc))/(np.amax(psc) - np.amin(psc))
#pims.append(p)

titles.append('outliers1')
psc = remove_outliers(asc,radius=5,threshold=100)
p = exposure.rescale_intensity(psc)
pims.append(p)

titles.append('outliers2')
psc = remove_outliers(asc,radius=3,threshold=20)
p = exposure.rescale_intensity(psc)
pims.append(p)

titles.append('top hat')
rob = white_tophat(asc,disk(2))
psc = asc - rob
p = exposure.rescale_intensity(psc)
pims.append(p)

titles.append('small objects')
bsc = asc > 0
rob = remove_small_objects(bsc,min_size=200,connectivity=20)
psc = asc * rob.astype('float')
p = exposure.rescale_intensity(psc)
pims.append(p)


#titles.append('mask')
#pims.append(rob.astype('float'))


#======================================================================
# plot

fig = plt.figure(figsize=[24,12])

ax = fig.add_subplot(2,3,1,projection=amap)
p = exposure.equalize_adapthist(asc)
pmap = sunpy.map.Map(p, x.header[idx])
pmap.plot(vmin = 0.0, vmax = 1.0)

frame = 2
for psc in pims:
   p = exposure.equalize_adapthist(psc)
   pmap = sunpy.map.Map(p, x.header[idx])
   ax = fig.add_subplot(2,3,frame,projection=pmap)
   pplot = pmap.plot(vmin = 0.0, vmax = 1.0, title=titles[frame-2])
   frame += 1

#======================================================================
plt.show()