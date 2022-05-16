"""
CIMP Background module
"""

# for warning / error statements; print red, yellow text to terminal
red = '\033[91m'
yellow = '\033[93m'
cend = '\033[0m'

class background:
    """
    The background class is used to compute a background from a series of     coronagraph images.  Typically the series will span 30 days of images     but it will work for other intervals.  The general procedure followed     for computing a background image is:
    1. Compute a daily median image for all images available for a given     day
    2. From the daily median images, construct a background image by     select the minimum positive (non-zero) value for each pixel.
    """

    def __init__(self, dir, N = 30):
        """
        dir is the directory containing the monthly data.  For now it is assumed that this directory contains subdirectories numbering from '01' to N, where N is the number of days, typically '30'.
        """

        self.dir = dir
        self.N = N
