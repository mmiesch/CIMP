
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import sunpy.visualization.colormaps as cm

from sunpy.io import fits

#dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/'
##file = dir+'03/32295364.fts'
#file = dir+'02/daily_median.fts'

#dir = '/home/mark.miesch/data/anny/CME0/pos-30/dcmer_030W_bang_0000_fits/tB/'
#file = dir+'frame_0050.fits'

#dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/'
dir = '/home/mark.miesch/data/lasco_monthly/c3/2014_01/'
file = dir+'background.fts'

#cmap_lasco_c2 = plt.get_cmap('soholasco2')
cmap_stereo_cor2 = plt.get_cmap('stereocor2')

data, header = fits.read(file)[0]

# experimenting with the model field

#header['CUNIT1'] = 'arcsec'
#header['CUNIT2'] = 'arcsec'
#
#header['CDELT1'] = 30/float(nx)
#header['CDELT2'] = 30/float(ny)
#
#header['CRPIX1'] = 242.5
#header['CRPIX2'] = 242.5

#print(f"{header['CRPIX1']} {header['CRPIX2']}")

dmap = sunpy.map.Map(data, header)
dmap.plot(cmap = cmap_stereo_cor2)

print(f" minmax: {np.nanmin(data)} {np.nanmax(data)}")

#plt.imshow(data, cmap = cmap_stereo_cor2)

plt.show()
