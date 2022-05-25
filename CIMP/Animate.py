"""
This module is for making movies!
"""

import os
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import sunpy.visualization.colormaps as cm

from astropy.time.core import Time
from CIMP import Snapshot as snap
from sunpy.io import fits

# for warning / error statements; print red, yellow text to terminal
red = '\033[91m'
yellow = '\033[93m'
cend = '\033[0m'

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

    def process(self, snap, background = 'ratio', method = 'none', \
                rmax = None):
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
            snap.subtract_background(rescale = False)
        else:
            snap.background_ratio(rescale = False)

        if method == 'enhance_mgn':
            snap.enhance(detail = 'mgn')
        elif method == 'enhance_fnrgf':
            snap.enhance(detail = 'fnrgf')

        if rmax is not None:
            snap.mask_annulus(rmax = rmax)

    def daymovie(self, background = 'ratio', method = 'None', \
                 scale = (0.0, 1.0), rmax = None, title = None, \
                 framedir = None, tolerance = None, diff_ratio = 10.0, \
                 resample = None, day = None):
        """
        loop over all valid files in a directory
        If you want to do nearest-neighbor interpolation on a regular grid, 
        then set `resample` equal to the integer number of desired grid points
        and `day` equal to the observation day as a string in iso time format, 
        e.g. "2022-05-25".  Frames with time stamps outside this day will be 
        ignored.
        """

        maps = []
        times = np.array([], dtype = 'float64')

        # get data from all valid files
        ref = None
        for file in os.listdir(self.dir):
            fpath = self.dir+'/'+file
            try:
                assert(os.path.isfile(fpath))
                assert("median" not in file)
                assert(fpath != self.bgfile)
                x = snap.snapshot(file = fpath, \
                    bgfile = self.bgfile, \
                    instrument = self.instrument, \
                    detector = self.detector)
                assert(x.nx == self.nx)
                assert(x.ny == self.ny)
                self.process(x, background = background, method = method, \
                             rmax = rmax)

                if x.valid(ref, tolerance, diff_ratio):
                    maps.append(x.map())
                    times = np.append(times, x.time.gps)
                    ref = x.data

            except:
                pass

        print(yellow+f"Nfiles = {len(maps)} {len(times)}"+cend)

        # optionally resample on to a regular time grid
        # set resample to be an integer equal to the number
        # of desired points in the regular grid
        if resample is not None:

            if day is None:

                # use the median as a robust estimator of the reference time
                # in case there are any invalid time stamps
                tref = np.nanmedian(times)
                times = times - tref

                # if the day is not specified, get a robust estimate of the 
                # time range by excluding any time stamps more than 20 hours 
                # away from the reference time
                vrange = 20.*3600.
                valid_times = np.ma.masked_where(np.abs(times) > vrange, times)
                tmin, tmax = (valid_times.min(), valid_times.max())

            else:
                tref = Time(day).gps
                tmin = 0.0
                tmax = Time(day+"T23:59:59.999").gps - tref
                times = times - tref
                tmin = np.max([tmin, np.nanmin(times)])
                tmax = np.min([tmax, np.nanmax(times)])

            # now create a regular array over the desired time range
            tgrid, dt = np.linspace(tmin, tmax, num = resample, endpoint = True, \
                                    retstep = True)

            # find nearest neighbor for each frame
            newmaps = []
            for t in tgrid:
                idx = (np.abs(times-t)).argmin()
                #print(f"{t} {times[idx]} {idx}")
                newmaps.append(maps[idx])
            maps = newmaps

        # make movie
        fig = plt.figure()
        frames = []
        for map in maps:
            if title is None:
                im = map.plot(cmap = self.cmap, vmin = scale[0], \
                              vmax = scale[1])
            else:
                im = map.plot(cmap = self.cmap, vmin = scale[0], \
                              vmax = scale[1], title = title)
            frames.append([im])
            if framedir is not None:
                frame = str(len(frames)).zfill(3)
                plt.savefig(framedir+f"/frame_{frame}.png")

        mov = animation.ArtistAnimation(fig, frames, interval = 50, blit = True,
              repeat = True, repeat_delay = 1000)

        print(yellow+f"Nframes = {len(frames)}"+cend)

        mov.save(self.outfile)

