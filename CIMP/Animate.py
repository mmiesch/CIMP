"""
This module is for making movies!
"""
import skvideo.io
import sunpy.visualization.colormaps as cm

from CIMP import Snapshot as snap

class movie:
    """
    Class for making movies from snapshots
    """

    def __init__(self, dir, bgfile = None, outfile = 'movie', \
                 instrument = None, detector = None, cmap = None):

        # parent directory containing image files
        self.dir = dir
        self.bgfile = bgfile
        self.outfile = outfile
        self.instrument = instrument
        self.detector = detector

        if cmap is None:
            self.cmap = plt.get_cmap('soholasco2')
        else:
            self.cmap = cmap

