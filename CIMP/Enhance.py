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

def bright_point_filter(im, threshold = 2.0, radius = 20):
    """
    Simpler version for just bright points
    """
    amed = median(im, disk(radius))
    rob = im > threshold * amed
    return np.where(rob, amed, im)

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

def rescale(im):
    return exposure.rescale_intensity(im, out_range=(0,1))

def equalize(im):
    ceq = exposure.equalize_adapthist(im)
    return rescale(ceq)

def detail(im, header, filter = 'mgn', instrument = None, detector = None):
    """
    specification of instrument and detector are only needed if the filter is nrgf or fnrgf.  It is used to determine the field of view.
    """

    if filter == 'mgn':
        """
        Multiscale Gaussian Noise filter (Morgan & Druckmuller 2014)
        """
        b = sunkit_image.enhance.mgn(im, h = 0.7, gamma = 1.5)

    elif filter == 'fnrgf':
        """
        Fourier Normaling Radial Gradient Filter (Druckmullerova et al 2011)
        """
        amap = sunpy.map.Map(im, header)
        bmap = fnrgf(amap, instrument, detector, order = 20, rmix = rmix)
        b = bmap.data

    elif filter == 'contrast':
        asc = rescale(a)
        b = enhance_contrast_percentile(asc, disk(2), p0=.1, p1=.9)
    else:
        b = a

    return rescale(b)

def denoise(im, filter = 'bregman'):

    if filter == 'tv':
        c = denoise_tv_chambolle(im, weight = 0.2)
    elif filter == 'bregman':
        c = denoise_tv_bregman(im)
    elif noise_removal == 'median':
        c = median(im, disk(1))
    elif filter == 'nl_means"':
        c = denoise_nl_means(a, patch_size = 4)
    else:
        print(yellow+"Warning: no noise filter applied")
        c = im

    return rescale(c)


