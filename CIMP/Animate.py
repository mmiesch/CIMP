"""
This module is for making movies!
"""

import os
import matplotlib.pyplot as plt
import sunpy.visualization.colormaps as cm

from CIMP import Snapshot as snap
from matplotlib.animation import FFMpegWriter
from sunpy.io import fits

class movie:
    """
    Class for making movies from snapshots
    """

    def __init__(self, dir, bgfile = None, outfile = 'movie.mp4', \
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

        metadata = dict(title="CIMP movie", 
                        instrument = self.instrument,
                        detector = self.detector)
        self.writer = FFMpegWriter(fps = 15, metadata = metadata)

        if bgfile is None:
            self.background = None

            # determine correct resolution from first fits file in director
            reffile = next(filter(lambda file: file[-5:] == ".fits" or file[-4:] == ".fts", flist2), None)

            data, header = fits.read(reffile)[0]
            self.nx = header['NAXIS1']
            self.ny = header['NAXIS2']

        else:
            data, header = fits.read(bgfile)[0]
            self.background = data
            self.nx = header['NAXIS1']
            self.ny = header['NAXIS2']

    def daymovie(self):
        """
        loop over all valid files in a directory
        """

        fig = plt.figure()

        # get data from all valid files
        self.nframes = 0
        for file in os.listdir(self.dir):
            fpath = self.dir+'/'+file
            try:
                assert("median" not in file)
                assert(fpath != self.bgfile)
                data, header = fits.read(fpath)[0]
                assert(header['NAXIS1'] == self.nx)
                assert(header['NAXIS2'] == self.ny)

                im = plt.imshow(data, cmap = self.cmap)
                frame = im.get_array()
                self.writer.writeFrame(frame)
                self.nframes += 1
            except:
                pass

        print(f"Nframes= {self.nframes}")

