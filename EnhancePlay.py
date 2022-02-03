"""
Use this to 
"""

from http.server import BaseHTTPRequestHandler
import numpy as np
from CIMP import Event as ev
from sunpy.net import attrs as a
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
import sunkit_image.enhance

import astroscrappy
import noisegate as ng

from sunkit_image.radial import fnrgf

from skimage import exposure
from skimage.filters import median
from skimage.filters.rank import enhance_contrast, enhance_contrast_percentile
from skimage.morphology import disk, remove_small_objects, remove_small_holes
from skimage.restoration import (denoise_tv_chambolle, denoise_bilateral,
                                 denoise_wavelet, estimate_sigma, denoise_nl_means)

from savitzkygolay import filter2D

import os
os.environ.pop("QT_QPA_PLATFORM_PLUGIN_PATH")

#---------------------------------------------------

tcase = 1

if tcase == 1:
    testcase = 1
    plotframe = 9
    scale = (0.0, 1000.0)

elif tcase == 2:
    testcase = 8
    plotframe = 3
    scale = (0.0, 1000.0)

else:
    print("specify a valid tcase")
    exit()    

#---------------------------------------------------
def remove_outliers(im, radius = 2, threshold = 50):
    medim = median(im,disk(radius))
    outliers = ( (im > medim + threshold) |
                 (im < medim - threshold) )
    out = np.where(outliers, medim, im)
    return out

def point_filter(im, threshold=None, min_size=64, connectivity=2):

    amag = np.absolute(im)

    if threshold is None:
        threshold = 0.1 * np.max(amag)

    bmag = amag > threshold
    rob = remove_small_objects(bmag,min_size=min_size,connectivity=connectivity)
    p = im * rob.astype('float')

    return p

def point_filter_med(im, threshold=2.0, radius=20):

    amag = np.absolute(im)
    amed = median(amag, disk(radius))
    rob = amag < (1.0+threshold) * amed
    p = im * rob.astype('float')
    return p

#---------------------------------------------------
# Define base image

x = ev.event.testcase(testcase)

a = x[plotframe]
amap = x.map(plotframe)

header = x.header[plotframe]

print(80*'-')
print(x)
print(80*'-')
print(repr(x))
print(80*'-')

print(f"Data range: {amap.min()} {amap.max()}")

#---------------------------------------------------
# Choose your battle

comp = (0,15)

tag = None

images = []
scales = []
titles = []

if comp.count(0) > 0:
    titles.append("Base image")
    images.append(a)
    scales.append(scale)

if comp.count(1) > 0:
   titles.append("NRGF")
   x.nrgf()
   b = x[plotframe]
   images.append(b)
   scales.append((0.0,4.0))

if comp.count(2) > 0:
    # current enhance() method in Event class
    titles.append("Event.enhance()")
    x.enhance(clip = scale)
    b = x[plotframe]
    images.append(b)
    scales.append(None)

if comp.count(3) > 0:
    # median filter with adaptive equalization
    titles.append("3 = clip + median + ad-eq-hist") 
  
    b = a.clip(min = scale[0], max = scale[1])
    p = exposure.rescale_intensity(b)
    b = median(p, disk(1))
    beq = exposure.equalize_adapthist(b)
    images.append(beq)
    scales.append((0.0,1.0))

if comp.count(4) > 0:
    # median filter with adaptive equalization
    titles.append("4 = 3 + en_con_per") 
  
    b = a.clip(min = scale[0], max = scale[1])
    p = exposure.rescale_intensity(b)
    b = median(p, disk(1))
    p = enhance_contrast_percentile(b, disk(2), p0=.1, p1=.9)
    b = exposure.rescale_intensity(p)
    beq = exposure.equalize_adapthist(b)
    images.append(beq)
    scales.append((0.0,1.0))

if comp.count(5) > 0:
    # noise gate filter with adaptive equalization
    titles.append("5 = ng + clip + ad-eq-hist") 
  
    N = 12
    dcube = np.zeros((N,a.shape[0],a.shape[1]))

    for i in np.arange(N):
        dcube[i,:,:] = x[i+1]

    print(f"dcube shape: {dcube.shape}")

    ng.noise_gate_batch(dcube, cubesize=(6,12,12), model='hybrid', factor = 2.0)
 
    b = a
    b[:,:] = dcube[plotframe,:,:]
    b = b.clip(min=scale[0], max=scale[1])

    p = exposure.rescale_intensity(b)
    beq = exposure.equalize_adapthist(p)
    images.append(beq)
    scales.append((0.0,1.0))

if comp.count(6) > 0:
    # noise gate filter with adaptive equalization
    titles.append("5 = ng + clip + median + ad-eq-hist") 
  
    N = 12
    dcube = np.zeros((N,a.shape[0],a.shape[1]))

    for i in np.arange(N):
        dcube[i,:,:] = x[i+1]

    print(f"dcube shape: {dcube.shape}")

    ng.noise_gate_batch(dcube, cubesize=(6,12,12), model='hybrid', factor = 2.0)
 
    b = a
    b[:,:] = dcube[plotframe,:,:]
    p = median(b, disk(1))
    b = p.clip(min=scale[0], max=scale[1])

    p = exposure.rescale_intensity(b)
    beq = exposure.equalize_adapthist(p)
    images.append(beq)
    scales.append((0.0,1.0))

if comp.count(7) > 0:
    # median filter with adaptive equalization
    titles.append("sunkit_image.enhance.mgn") 
  
    b = a.clip(min = scale[0], max = scale[1])
    
    p = sunkit_image.enhance.mgn(b); tag='a'
    #p = sunkit_image.enhance.mgn(b,h=0.1); tag='b'
    #p = sunkit_image.enhance.mgn(b,h=1.0); tag='c'
    #p = sunkit_image.enhance.mgn(b,h=1.0,k=2.0); tag='d'
    #p = sunkit_image.enhance.mgn(b,h=1.0,k=0.1); tag='e'
    #p = sunkit_image.enhance.mgn(b,h=0.8,gamma=2.5); tag='f'
    p = sunkit_image.enhance.mgn(b,k=0.8,h=0.7,sigma=[1.25, 2.5, 5, 10, 20, 40],
                                weights=[0.5,0.5,0.5,0.5,0.5,10.0],gamma=2); tag='g'

    b = exposure.rescale_intensity(p)
    #beq = exposure.equalize_adapthist(b)
    beq = exposure.equalize_hist(b)
    images.append(beq)
    scales.append(None)

if comp.count(8) > 0:
    titles.append("astroscrappy") 
  
    b = a.clip(min = scale[0], max = scale[1])
    p = exposure.rescale_intensity(b)

    mask, psc = astroscrappy.detect_cosmics(p, sigclip=2, objlim=2, readnoise=4, verbose=True)

    b = exposure.rescale_intensity(psc)
    beq = exposure.equalize_adapthist(b)
    images.append(beq)
    scales.append(None)

if comp.count(9) > 0:
    titles.append("outliers") 
  
    b = a.clip(min = scale[0], max = scale[1])

    p = exposure.rescale_intensity(b)
    psc = remove_outliers(p,radius=3,threshold=20)
    p = exposure.rescale_intensity(psc)

    beq = exposure.equalize_adapthist(p)
    images.append(beq)
    scales.append([0,1])

if comp.count(10) > 0:
    titles.append("remove small objects") 
  
    mx = scale[1]

    #b = point_filter(a, threshold = 0.1*mx, min_size=64, connectivity=2); tag='a'
    #b = point_filter(a, threshold = 0.1*mx, min_size=64, connectivity=20); tag='b'
    #b = point_filter(a, threshold = 0.1*mx, min_size=10, connectivity=20); tag='c'
    b = point_filter(a, threshold = 0.05*mx, min_size=64, connectivity=20); tag='c'

    images.append(b)
    scales.append(scale)

if comp.count(11) > 0:
    titles.append("Savitzky-Golay") 
  
    b = a.clip(min = scale[0], max = scale[1])

    p = filter2D(b,mu=5,poly=3,order=1)

    b = exposure.rescale_intensity(p)
    #beq = exposure.equalize_adapthist(p)
    images.append(b)
    scales.append([.2,.7])

if comp.count(12) > 0:
    titles.append("FNRGF") 
    fov = {
       'lasco-c2'   : (1.5,  6.0),
       'lasco-c3'   : (3.7, 30.0),
       'secchi-cor1': (1.5,  4.0),
       'secchi-cor2': (3.0, 15.0)
    }
    
    b = a.clip(min = scale[0], max = scale[1])
    bmap = sunpy.map.Map(b,header)

    myfov = fov[x.instrument.lower()+'-'+x.detector.lower()]
        
    edges = equally_spaced_bins(myfov[0], myfov[1])
    edges *= u.R_sun

    map = radial.nrgf(bmap, edges)

    p = map.data

    b = exposure.rescale_intensity(p)
    #beq = exposure.equalize_adapthist(p)
    images.append(p)
    scales.append((0,1))

if comp.count(13) > 0:
    titles.append("remove small objects + clip + adeq") 
  
    mx = scale[1]

    b = point_filter(a, threshold = 0.05*mx)

    bclip = b.clip(min = scale[0], max = scale[1])
    p = exposure.rescale_intensity(bclip)
    beq = exposure.equalize_adapthist(p)
    images.append(beq)
    scales.append([0,1])

if comp.count(14) > 0:
    titles.append("point_filter_med") 
  
    b = point_filter_med(a); tag='a'
    
    images.append(b)
    scales.append(scale)

if comp.count(15) > 0:
    titles.append("point_filter_med + clip + adeq") 
  
    b = point_filter_med(a)
    
    p = exposure.rescale_intensity(b)
    beq = exposure.equalize_adapthist(p)

    images.append(beq)
    scales.append([0,1])
       

# ===================

fig = plt.figure(figsize=[16,8])

map1 = sunpy.map.Map(images[0],header)
print(f"image 1 range: {map1.min()} {map1.max()}")
ax = fig.add_subplot(1,2,1,projection=map1)
if scales[0] is None:
    map1.plot(title=titles[0])
else:
   print(f"image 1 scale: {scales[0][0]} to {scales[0][1]}")
   map1.plot(vmin=scales[0][0], vmax=scales[0][1], title=titles[0])

map2 = sunpy.map.Map(images[1],header)
print(f"image 2 range: {map2.min()} {map2.max()}")
ax = fig.add_subplot(1,2,2,projection=map2)
if scales[1] is None:
    map2.plot(title=titles[1])
else:
   print(f"image 2 scale: {scales[1][0]} to {scales[1][1]}")
   map2.plot(vmin=scales[1][0], vmax=scales[1][1], title=titles[1])

# ===================
# save to a file
# ===================
dir = '/home/mark.miesch/Products/image_processing/CIMP/images/Enhance_comp/'

if tag is None:
   file = dir + f"Enhance_t{testcase}_{comp[0]}_vs_{comp[1]}.png"
else:
   file = dir + f"Enhance_t{testcase}_{comp[0]}_vs_{comp[1]}_{tag}.png"

plt.savefig(file)

plt.show()