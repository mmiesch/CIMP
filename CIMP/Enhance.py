"""
CIMP Enhance module
"""

import astropy.units as u
import sunkit_image.radial as radial
import sunpy.map

from sunkit_image.utils import equally_spaced_bins

fov = {
    'lasco-c2'   : (1.5,  6.0),
    'lasco-c3'   : (3.7, 30.0),
    'secchi-cor1': (1.5,  4.0),
    'secchi-cor2': (3.0, 15.0)
}

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
