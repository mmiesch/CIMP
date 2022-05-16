
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import sunpy.io
import sunpy.visualization.colormaps as cm


dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/'
file = dir+'03/32295364.fts'

#dir = '/home/mark.miesch/data/anny/CME0/pos-30/dcmer_030W_bang_0000_fits/tB/'
#file = dir+'frame_0050.fits'

#cmap_lasco_c2 = plt.get_cmap('soholasco2')
cmap_stereo_cor2 = plt.get_cmap('stereocor2')

data, header = sunpy.io.fits.read(file)[0]

nx = data.shape[0]
ny = data.shape[0]

# simple r^2 scaling
x0 = 0.5 * nx
y0 = 0.5 * ny

for i in np.arange(0,nx):
    for j in np.arange(0,ny):
        r2 = (float(i)-x0)**2 + (float(j)-y0)**2
        if r2 > 1:
            data[i,j] = data[i,j] * r2

#for h in header.keys():
#    print(f"{h} : {header[h]}")

#try:
#    dmap = sunpy.map.Map(data, header)
#    dmap.plot()
#except:
#    plt.imshow(data)

header['CUNIT1'] = 'arcsec'
header['CUNIT2'] = 'arcsec'

header['CDELT1'] = 30/float(nx)
header['CDELT2'] = 30/float(ny)

header['CRPIX1'] = 242.5
header['CRPIX2'] = 242.5

print(f"{header['CRPIX1']} {header['CRPIX2']}")

dmap = sunpy.map.Map(data, header)
dmap.plot(cmap = cmap_stereo_cor2)

#plt.imshow(data, cmap = cmap_stereo_cor2)

#plt.show()