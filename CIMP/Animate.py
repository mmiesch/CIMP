"""
This module is for making movies!
"""

import os
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import sunpy.visualization.colormaps as cm

from astropy.io.fits import ImageDataDiff
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
                rmin = 0.0, rmax = None):
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

        if method == 'enhance_mgn':
            snap.enhance(point = 'omr', detail = 'mgn', noise_filter = 'omr')
        elif method == 'enhance_fnrgf':
            snap.enhance(point = 'omr', detail = 'fnrgf', noise_filter = 'omr')

        snap.mask_annulus(rmin = rmin, rmax = rmax)

    def valid(self, image, ref, tolerance = None, diff_ratio = 100.0, \
              file = ''):
        """
        This identifies and flags corrupted images by comparing the data to a reference image that is passed as ref (a 2D numpy array).  In a movie sequence, ref would typically be the previous frame.
        Images are flagged as corrupted if more than `diff_ratio` percent of their pixels are different from the reference image, within a specified `tolerance`.  The default diff_ratio is set to be 100%,
        which means that, by default, no frames are declared invalid.
        The default tolerance is currently set to be 0.5 times the median non-zero value of the reference image.
        """

        assert(image.shape == ref.shape)

        if tolerance is None:
            valid_pix = np.ma.masked_equal(ref, 0.0, copy = False)
            tolerance = 0.5*np.nanmedian(np.absolute(valid_pix))

        d = ImageDataDiff(image, ref, rtol = tolerance)

        if (100*d.diff_ratio > diff_ratio):
            print(yellow+f"Corrupted image {file} {tolerance} {diff_ratio} {d.diff_ratio*100} {np.min(ref)} {np.max(ref)}" + cend)
            return False
        else:
            return True

    def daymovie(self, background = 'ratio', method = 'None', \
                 scale = (0.0, 1.0), rmin = None, rmax = None, title = None, \
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

        # first pass: get data from all valid files
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
                             rmin = rmin, rmax = rmax)
                maps.append(x.map())
                times = np.append(times, x.time.gps)

            except:
                print(red+f"Skipping {file}"+cend)
                pass

        Nmaps = len(maps)
        print(yellow+f"Nfiles = {Nmaps}"+cend)

        # second pass: remove corrupted images
        valid_maps = []
        valid_times = []
        for idx in np.arange(Nmaps):
            i1 = np.max([0, idx-2])
            i2 = np.min([idx+3, Nmaps])
            Nref = i2 - i1
            refimages = np.empty((self.nx, self.ny, Nref), dtype = 'float32')
            for ridx in np.arange(Nref):
                refimages[:,:,ridx] = maps[ridx+i1].data
            ref = np.nanmedian(refimages,axis=2)
            if self.valid(maps[idx].data, ref, tolerance, diff_ratio, \
                          file = maps[idx].fits_header['FILENAME']):
                valid_maps.append(maps[idx])
                valid_times.append(times[idx])

        # free up some memory before resampling
        maps = valid_maps
        times = valid_times
        refimages = 0
        valid_maps = 0
        valid_times = 0
        N_valid_files = len(maps)

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
                print(yellow+f"{tmin} {t} {times[idx]} {tmax} {idx}"+cend)
                newmaps.append(maps[idx])
            maps = newmaps

        # make movie
        fig = plt.figure()
        frames = []
        for map in maps:
            im = plt.imshow(map.data, cmap=self.cmap, vmin = scale[0], \
                            vmax = scale[1], origin='lower')
            if title is not None:
                plt.title(title)
            frames.append([im])
            if framedir is not None:
                frame = str(len(frames)).zfill(3)
                plt.savefig(framedir+f"/frame_{frame}.png")

        mov = animation.ArtistAnimation(fig, frames, interval = 50, blit = True,
              repeat = True, repeat_delay = 1000)

        print(yellow+f"Number of valid files = {N_valid_files}"+cend)
        print(yellow+f"Number of frames = {len(frames)}"+cend)

        mov.save(self.outfile)

