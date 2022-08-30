
import os

from astropy.io import fits

#dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04/10'
dir = '/home/mark.miesch/data/lasco_monthly/c3/2014_01/11'

print(dir)

for file in os.listdir(dir):
    fpath = dir+'/'+file
    hdu = fits.open(fpath)
    print(hdu[0].header['TIME-OBS'])
    hdu.close()


