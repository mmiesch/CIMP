"""
This is a driver routine to test the computation and subtraction of the background.
"""

import matplotlib.pyplot as plt
import numpy as np

from CIMP import Background as bg

print("Starting")

# for snapshot testcase 1
dir = '/home/mark.miesch/data/lasco_monthly/c3/2012_04'

# for snapshot testcase 2
#dir = "/home/mark.miesch/data/lasco_monthly/c3/2014_01"

# another lasco test set
#dir = "/home/mark.miesch/data/lasco_monthly/c3/2021_05"

# for full-res STEREO-A L1 data
#dir = "/home/mark.miesch/sunpy/data/secchi_cor2/L1/2012/09"

normalize = False
lasco_correction = True

bg = bg.background(dir)

bg.daily_medians(normalize = normalize, lasco_correction = lasco_correction)
bg.minimize_medians()
bg.write_background()


