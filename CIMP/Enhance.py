"""
CIMP Enhance module
"""

import astropy.units as u
import numpy as np
import sunkit_image.enhance
import sunkit_image.radial as radial
import sunpy.map

from skimage import exposure
from skimage.filters import median
from skimage.filters.rank import enhance_contrast_percentile
from skimage.morphology import disk
from skimage.restoration import (denoise_tv_chambolle, denoise_tv_bregman,
                                 denoise_nl_means)
from sunkit_image.utils import equally_spaced_bins

fov = {
    'lasco-c2'   : (1.5,  6.0),
    'lasco-c3'   : (3.7, 30.0),
    'secchi-cor1': (1.5,  4.0),
    'secchi-cor2': (3.0, 15.0),
    'cme-model1' : (6.0, 120.0)
}

def point_filter(im, threshold = 2.0, radius = 20):
    amag = np.absolute(im)
    amag -= np.min(amag)
    amed = median(amag, disk(radius))
    rob = amag < (1.0+threshold) * amed
    p = im * rob.astype('float')
    return np.where(amag > 0, p, 0.0)

def nrgf(imap, instrument = 'lasco', detector = 'c3'):

    myfov = fov[instrument.lower()+'-'+detector.lower()]

    edges = equally_spaced_bins(myfov[0], myfov[1])
    edges *= u.R_sun

    return radial.nrgf(imap, edges)

def fnrgf(imap, instrument = 'lasco', detector = 'c3', order = 20, rmix = [1,15]):
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

def enhance(imap, instrument = 'lasco', detector = 'c3', clip = None,
            noise_filter = 'bregman', detail = 'mgn', rmix = [1,15]):
    """
    A combination of actions that works pretty well
    """

    # filter out bright points
    a = point_filter(imap.data)

     # contrast stretching via clipping
    if clip is not None:
        a = a.clip(min = clip[0], max = clip[1])

    # various techniques to bring out detail
    if detail == 'mgn':
        """
        Multiscale Gaussian Noise filter (Morgan & Druckmuller 2014)
        """
        b = sunkit_image.enhance.mgn(a, h = 1.0, gamma = 1.5)

    elif detail == 'fnrgf':
        """
        Fourier Normaling Radial Gradient Filter (Druckmullerova et al 2011)
        """
        amap = sunpy.map.Map(a, imap.meta)
        bmap = fnrgf(amap, instrument, detector, order = 20, rmix = rmix)
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

    return c

