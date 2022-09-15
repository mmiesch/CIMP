"""
check the lasco exposure time correction
do this before running:
export NRL_LIB=/home/mark.miesch/data/lasco_cal
"""

import os
import sys
from astropy.io import fits

sys.path.append('/home/mark.miesch/Products/image_processing/NRL_python')
from lasco_utils import get_exp_factor

caldir = '/home/mark.miesch/data/lasco_cal/idl/expfac/data'

#infile = '/home/mark.miesch/data/lasco_monthly/c3/2014_01/20/33385905.fts'

#infile = '/home/mark.miesch/data/lasco_monthly/c3/2021_05/16/33684538.fts'

#dir = '/home/mark.miesch/data/lasco_monthly/c3/2014_01/20'

#dir='/home/mark.miesch/Products/image_processing/ATBD/data/lasco_c3/L2proxy_2014_01'
#dir='/home/mark.miesch/Products/image_processing/ATBD/data/lasco_c3/L2proxy_2021_05'
#dir='/home/mark.miesch/Products/image_processing/ATBD/data/lasco_c3/L2proxy_2012_04'

dir='/home/mark.miesch/data/lasco_monthly/c3/2021_05/11'

#------------------------------------------------------------------------------

for d in os.listdir(dir):
   infile = dir+'/'+d
   hdu = fits.open(infile)
   fac, bias = get_exp_factor(hdu[0].header, dir0 = caldir)
   print(f"{fac} {bias}")

