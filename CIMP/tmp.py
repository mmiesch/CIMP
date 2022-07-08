
import astropy.units as u
import sunkit_image.radial as radial
from sunpy.map import Map
from sunkit_image.utils import equally_spaced_bins

rmin = 3.0
rmax = 15.0
order = 20
Bmap = Map(image, header)
edges = equally_spaced_bins(rmin, rmax)
edges *= u.R_sun
coefs = radial.set_attenuation_coefficients(order)
Rmap = radial.fnrgf(Bmap, edges, order, coefs, ratio_mix = [1,15])
result = Rmap.data

