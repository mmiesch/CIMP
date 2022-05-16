
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import sunpy.visualization.colormaps as cm

from sunpy.io import fits


dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/'

cmap_stereo_cor2 = plt.get_cmap('stereocor2')

Ndays = 30

for d in np.arange(1,Ndays+1):
    day = f'{d:02d}'
    file = dir + day + '/daily_median.fts'
    data, header = fits.read(file)[0]
    dmap = sunpy.map.Map(data, header)
    dmap.plot(cmap = cmap_stereo_cor2)
    plt.savefig('/home/mark.miesch/tmp/'+day+".png",format='png')

