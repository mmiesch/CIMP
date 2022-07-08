


import astropy.units as u
import sunkit_image.radial as radial
from skimage import exposure
from sunpy.map import Map
from sunkit_image.utils import equally_spaced_bins

kmax = 20
scaled_image = exposure.rescale_intensity(image.clip(min = 1.0, max = 1.3), out_range=(0,1))
Bmap = Map(scaled_image, header)
edges = equally_spaced_bins(2.45, 15.0) * u.R_sun
coefs = radial.set_attenuation_coefficients(kmax)
Rmap = radial.fnrgf(Bmap, edges, kmax, coefs, ratio_mix = [4,1])
result = Rmap.data



