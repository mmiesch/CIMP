"""
check the lasco exposure time correction
"""

import sys

from astropy.io import fits

sys.path.append('/home/mark.miesch/Products/image_processing/NRL_python')
from lasco_utils import get_exp_factor

caldir = '/home/mark.miesch/data/lasco_cal'


