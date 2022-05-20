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

    Current options for 

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

    def process(self, snap, background = 'ratio', method = 'none'):
        """
        Image processing to be performed for each snapshot

        Current options for background are:
        'subtract': subtract background
        'ratio': ratio of image over background

        Current options for method are:
        'none': no processing
        'enhance mgn': enhance method with mgn
        'enhance fnrgf': enhance method with fnrgf
        """

        if background == 'subtract':
            snap.subtract_background()
        else:
            snap.background_ratio()

        if background == 'enhance mgn':
            snap.enhance(detail = 'mgn')
        elif background == 'enhance fnrgf':
            snap.enhance(detail = 'fnrgf')

    def daymovie(self, background = 'ratio', method = 'None', \
                 scale = (0.0, 1.0)):
        """
        loop over all valid files in a directory
        """

        fig = plt.figure()
        frames = []

        # get data from all valid files
        for file in os.listdir(self.dir):
            fpath = self.dir+'/'+file
            try:
                assert("median" not in file)
                assert(fpath != self.bgfile)
                x = snap.snapshot(file = fpath, \
                    bgfile = self.bgfile, \
                    instrument = self.instrument, \
                    detector = self.detector)
                assert(x.nx == self.nx)
                assert(x.ny == self.ny)
                self.process(x, background = background, method = method)
                im = x.map().plot(cmap = self.cmap, vmin = scale[0], \
                                  vmax = scale[1])
                frames.append([im])
            except:
                pass

        mov = animation.ArtistAnimation(fig, frames, interval = 50, blit = True,
              repeat = True, repeat_delay = 1000)

        print(f"Nframes= {len(frames)}")

        mov.save(self.outfile)


