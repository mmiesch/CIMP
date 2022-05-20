"""
This module is for making movies!
"""

import os
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import sunpy.visualization.colormaps as cm

from CIMP import Snapshot as snap
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
        frames = []

        # get data from all valid files
        for file in os.listdir(self.dir):
            fpath = self.dir+'/'+file
            print(f"{file} {self.nx} {self.ny}")
            #try:
            #assert("median" not in file)
            #assert(fpath != self.bgfile)
            data, header = fits.read(fpath)[0]
            #assert(header['NAXIS1'] == self.nx)
            #assert(header['NAXIS2'] == self.ny)
            im = plt.imshow(data, cmap = self.cmap)
            frames.append([im])
            #except:
            #    pass

        mov = animation.ArtistAnimation(fig, frames, interval = 50, blit = True,
              repeat = True, repeat_delay = 1000)

        print(f"Nframes= {len(frames)}")

        mov.save(self.outfile)


