"""
CIMP Enhance module
"""

import astropy.units as u
import numpy as np
import sunkit_image.enhance
import sunkit_image.radial as radial
import sunpy.map

from scipy.signal import convolve2d
from skimage import exposure
from skimage.filters import median
from skimage.filters.rank import enhance_contrast_percentile
from skimage.measure import block_reduce
from skimage.morphology import disk, opening, erosion, reconstruction
from skimage.restoration import (denoise_tv_chambolle, denoise_tv_bregman,
                                 denoise_nl_means, denoise_wavelet)
from sunkit_image.utils import equally_spaced_bins

# for warning / error statements; print red, yellow text to terminal
red = '\033[91m'
yellow = '\033[93m'
cend = '\033[0m'

fov = {
    'lasco-c2'   : (1.5,  6.0),
    'lasco-c3'   : (3.5, 30.0),
    'secchi-cor1': (1.5,  4.0),
    'secchi-cor2': (2.45, 15.0),
    'cme-model1' : (6.0, 120.0)
}

def rescale(im):
    return exposure.rescale_intensity(im, out_range=(0,1))

def point_filter(im, threshold = 2.0, radius = 20, rescaleim = True):
    amag = np.absolute(im)
    amag -= np.min(amag)
    amed = median(amag, disk(radius))
    rob = amag < (1.0+threshold) * amed
    p = im * rob.astype('float')
    a = np.where(amag > 0, p, 0.0)
    if rescaleim:
        return rescale(a)
    else:
        return a

def bright_point_filter(im, threshold = 2.0, radius = 20, rescaleim = False):
    """
    Simpler version for just bright points
    """
    amed = median(im, disk(radius))
    rob = im > threshold * amed
    a = np.where(rob, amed, im)
    if rescaleim:
        return rescale(a)
    else:
        return a

def morph_opening(im, radius = 3, rescaleim = True):
    """
    Point filter based on mathematical morphology
    """
    a = opening(im, disk(radius))
    if rescaleim:
        return rescale(a)
    else:
        return a

def omr(im, radius = 3, rescaleim = True):
    """
    OMR = Opening by Morphological Reconstruction
    This is similar to morphological opening but uses the original image as a mask to mitigate the apparent smoothing that comes with a standard opening procedure.
    """
    a = erosion(im, disk(radius))
    b = reconstruction(a, im, method = 'dilation', footprint = disk(radius))

    if rescaleim:
        return rescale(b)
    else:
        return b

def nrgf(imap, instrument = 'lasco', detector = 'c3'):

    myfov = fov[instrument.lower()+'-'+detector.lower()]

    edges = equally_spaced_bins(myfov[0], myfov[1])
    edges *= u.R_sun

    return radial.nrgf(imap, edges)

def fnrgf(imap, instrument = 'lasco', detector = 'c3', order = 20, rmix = [4,1]):
    """
    Fourier Normaling Radial Gradient Filter (Druckmullerova et al 2011)
    """
    myfov = fov[instrument.lower()+'-'+detector.lower()]

    edges = equally_spaced_bins(myfov[0], myfov[1])
    edges *= u.R_sun

    coefs = radial.set_attenuation_coefficients(order)

    return radial.fnrgf(imap, edges, order, coefs, ratio_mix = rmix)

def powerlaw(im, n = 2.0, center = None):
    """
    Offsetting the radial gradient with a simple power law
    """

    nover2 = n / 2.0

    nx = im.shape[0]
    ny = im.shape[1]

    norm = (0.25*float(nx))**2

    if center is None:
        center = (0.5*float(nx), 0.5*float(ny))

    f = np.zeros((nx,ny),dtype='float32')

    for i in np.arange(0,nx):
        for j in np.arange(0,ny):
            r2 = (float(i)-center[0])**2 + (float(j)-center[1])**2
            if r2 > 0:
                f[i,j] = im[i,j] * (r2/norm)**nover2

    return f

def clip(im, min = None, max = None, rescale_output = False):
    if min is None:
        min = np.min(im)
    if max is None:
        max = np.max(im)

    a = im.clip(min = min, max = max)

    if rescale_output:
        return rescale(a)
    else:
        return a

def equalize(im):
    ceq = exposure.equalize_adapthist(im)
    return rescale(ceq)

def detail(im, header, filter = 'mgn', instrument = None, detector = None, \
           params = None):
    """
    specification of instrument and detector are only needed if the filter is nrgf or fnrgf.  It is used to determine the field of view.
    """

    rescale_output = True

    if filter == 'mgn':
        """
        Multiscale Gaussian Noise filter (Morgan & Druckmuller 2014)
        """
        if params is None:
            h = 0.8
            gamma = 1.5
        else:
            h = params[0]
            gamma = params[1]

        b = sunkit_image.enhance.mgn(im, h = h, gamma = gamma)

    elif filter == 'nrgf':
        """
        Normaling Radial Gradient Filter (Morgan et al 2006)
        """
        amap = sunpy.map.Map(im, header)
        bmap = nrgf(amap, instrument, detector)
        b = bmap.data

    elif filter == 'fnrgf':
        """
        Fourier Normaling Radial Gradient Filter (Druckmullerova et al 2011)
        """
        amap = sunpy.map.Map(im, header)
        bmap = fnrgf(amap, instrument, detector)
        b = bmap.data

    elif filter == 'contrast':
        asc = rescale(a)
        b = enhance_contrast_percentile(asc, disk(2), p0=.1, p1=.9)

    else:
        print(yellow+"Warning: no detail enhancement applied"+cend)
        b = im
        rescale_output = False

    if rescale_output:
        return rescale(b)
    else:
        return b

def denoise(im, filter = 'bregman'):

    rescale_output = True

    if filter == 'tv':
        c = denoise_tv_chambolle(im, weight = 0.1)
    elif filter == 'bregman':
        c = denoise_tv_bregman(im, weight = 10)
    elif filter == 'median':
        c = median(im, disk(1))
    elif filter == 'nl_means"':
        c = denoise_nl_means(im, patch_size = 4)
    elif filter == 'omr':
        c = omr(im, rescaleim = False)
    elif filter == 'atrous':
        c = atrous(im)
    elif filter == 'wavelet':
        c = denoise_wavelet(im,mode='hard',wavelet='haar')
    else:
        print(yellow+"Warning: no noise filter applied"+cend)
        c = im
        rescale_output = False

    if rescale_output:
        return rescale(c)
    else:
        return c

def atrous(im):
    """
    a trous wavelet filtering
    """

    kernel_b3_spline = np.array([
        [1 / 256, 1 / 64, 3 / 128, 1 / 64, 1 / 256],
        [1 / 64, 1 / 16, 3 / 32, 1 / 16, 1 / 64],
        [3 / 128, 3 / 32, 9 / 64, 3 / 32, 3 / 128],
        [1 / 64, 1 / 16, 3 / 32, 1 / 16, 1 / 64],
        [1 / 256, 1 / 64, 3 / 128, 1 / 64, 1 / 256]
    ])

    c0 = im
    c1 = convolve2d(c0, kernel_b3_spline, mode='same', boundary='fill')

    return c1

def mask_annulus(im, rmin = 0.0, rmax = np.inf, missingval = 0.0):
    """
    This sets the pixels inside rmin and/or outside rmax to the missing value (default 0)
    """
    # annular radii in pixels (squared)
    nn = 0.5 * float(np.min((im.shape[0], im.shape[1])))
    rr1 = (rmin*nn)**2
    rr2 = (rmax*nn)**2

    # center of image
    x0 = 0.5*float(im.shape[0])
    y0 = 0.5*float(im.shape[1])

    x = np.tile(np.arange(im.shape[0]),(im.shape[1],1)).T
    y = np.tile(np.arange(im.shape[1]),(im.shape[0],1))
    rr = (x-x0)**2 + (y-y0)**2
    mask = (rr <= rr1) | (rr >= rr2)
    return np.where(mask, missingval, im)

def downsample(im):
    """
    This downsamples an image by a factor of two in each dimension.  So, for example, it will convert a 2048 x 2048 image into a 1024 x 1024 image.  The method it uses to downsample is to replace each 2 x 2 block in the source image by the median of those four values.  So, it serves as a simple point filter while downsampling.
    """

    return block_reduce(im, block_size = 2, func = np.nanmedian, cval = np.nanmin(im))
