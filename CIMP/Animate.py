"""
This module is for making movies!
"""

import os
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import noisegate as ng
import numpy as np
import sunpy.map
import sunpy.visualization.colormaps as cm

from astropy.io import fits
from astropy.io.fits import ImageDataDiff
from astropy.time.core import Time
from CIMP import Snapshot as snap
from CIMP import Enhance

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

            # determine correct resolution from first fits file in directory
            flist = os.listdir(self.dir)
            rfile = next(filter(lambda file: file[-5:] == ".fits" or file[-4:] == ".fts", flist), None)
            reffile = self.dir+'/'+rfile
            print(yellow+f"Reference file: {reffile}"+cend)

            hdu = fits.open(reffile)[0]
            self.nx = hdu.header['NAXIS1']
            self.ny = hdu.header['NAXIS2']

        else:
            hdu = fits.open(bgfile)[0]
            self.background = hdu.data
            self.nx = hdu.header['NAXIS1']
            self.ny = hdu.header['NAXIS2']

    def process(self, snap, background = 'ratio', method = 'none', \
                downsample = False, clip = None, rmin = 0.0, rmax = None, \
                params = None):
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
            snap.mask_background(rmin = rmin, rmax = rmax, nonzero = False)
            snap.subtract_background(rescale=False)
        else:
            snap.mask_background(rmin = rmin, rmax = rmax, nonzero = True)
            snap.background_ratio(rescale=False)

        snap.mask_annulus(rmin = rmin, rmax = rmax)

        if downsample:
            snap.downsample()

        if method == 'enhance_mgn':
            snap.enhance(clip = clip, point = 'omr', detail = 'mgn', noise_filter = 'none', params = params)
        elif method == 'enhance_fnrgf':
            snap.enhance(clip = clip, point = 'omr', detail = 'fnrgf', noise_filter = 'none')

        # hit it with another mask after processing
        snap.mask_annulus(rmin=rmin, rmax = rmax, missingval = np.nanmin(snap.data))

    def noise_gate(self, maps, cubesize = (3, 12, 12), model = 'hybrid',
                   factor = 2.0):
        """
        Noise Gate filter from DeForest, C.E. 2017, ApJ, 838:155 (10pp)
        """

        nframes = len(maps)

        print(f"Applying noise gate filter: {cubesize} {model} {factor}")
        dcube = np.zeros((nframes, maps[0].data.shape[0],
                                   maps[0].data.shape[1]))

        for i in np.arange(1, nframes):
            dcube[i-1,:,:] = maps[i].data

        ng.noise_gate_batch(dcube, cubesize=cubesize, model=model,
                            factor = factor)

        for i in np.arange(1, nframes):
            maps[i] = sunpy.map.Map(dcube[i-1,:,:], maps[i].fits_header)

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

    def daymovie(self, background = 'ratio', method = 'None', clip = None, \
                 scale = (0.0, 1.0), rmin = None, rmax = None, title = None, \
                 framedir = None, tolerance = None, diff_ratio = 10.0, \
                 resample = None, day = None, noisegate = False, \
                 downsample = False, params = None):
        """
        loop over all valid files in a directory
        If you want to do nearest-neighbor interpolation on a regular grid, 
        then set `resample` equal to the integer number of desired grid points
        and `day` equal to the observation day as a string in iso time format, 
        e.g. "2022-05-25".  Frames with time stamps outside this day will be 
        ignored.
        """

        maps = []
        files = []
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
                             downsample = downsample, clip = clip, \
                             rmin = rmin, rmax = rmax, params = params)
                maps.append(x.map())
                times = np.append(times, x.time.gps)
                files.append(file)
            except:
                print(red+f"Skipping {file}"+cend)
                pass

        Nmaps = len(maps)
        print(yellow+f"Nfiles = {Nmaps}"+cend)

        # second pass: remove corrupted images
        valid_maps = []
        valid_files = []
        valid_times = []
        nx = maps[0].data.shape[0]
        ny = maps[0].data.shape[1]
        for idx in np.arange(Nmaps):
            i1 = np.max([0, idx-2])
            i2 = np.min([idx+3, Nmaps])
            Nref = i2 - i1
            refimages = np.empty((nx, ny, Nref), dtype = 'float32')
            for ridx in np.arange(Nref):
                refimages[:,:,ridx] = maps[ridx+i1].data
            ref = np.nanmedian(refimages,axis=2)
            if self.valid(maps[idx].data, ref, tolerance, diff_ratio, \
                          file = maps[idx].fits_header['FILENAME']):
                valid_maps.append(maps[idx])
                valid_files.append(files[idx])
                valid_times.append(times[idx])

        # free up some memory before resampling
        maps = valid_maps
        files = valid_files
        times = valid_times
        refimages = 0
        valid_maps = 0
        valid_times = 0
        N_valid_files = len(maps)

        # show which file corresponds to each frame (before resampling)
        for idx in np.arange(0,len(files)):
            print(f"frame {idx} = {files[idx]}")

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

        # optionally apply a noise gate filter
        if noisegate:
            self.noise_gate(maps)

        # make movie
        fig = plt.figure()
        frames = []
        for map in maps:
            im = plt.figimage(map.data, cmap=self.cmap, vmin = scale[0], \
                            vmax = scale[1], origin='lower', resize=True)
            if title is not None:
                plt.title(title)
            frames.append([im])
            if framedir is not None:
                frame = str(len(frames)).zfill(3)
                plt.savefig(framedir+f"/frame_{frame}.png")

        mov = animation.ArtistAnimation(fig, frames, interval = 50, blit = False,
              repeat = True, repeat_delay = 1000)

        print(yellow+f"Number of valid files = {N_valid_files}"+cend)
        print(yellow+f"Number of frames = {len(frames)}"+cend)

        mov.save(self.outfile)


    def dirmovie(self, scale = (0.0, 1.0), title = None, \
                 denylist = None, tolerance = None, diff_ratio = 10.0, \
                 framedir = None, fitsdir = None):
        """
        This is a simplified version of daymovie.  Like daymovie, it loops over
        all fits files in a directory to make a movie.  However, it does not
        apply any processing and it does not interpolate onto a regular time grid.
        But, it does do some qc and, optionally, you can write the qc-filtered
        files to outdir.  This is useful as a pre-processing stage before applying
        a noisegate filter.  So, the goal is to produce a series of files with a
        similar noise profile (filtering out degraded images) that can be modeled
        by noisegate.
        """

        images = []
        headers = []

        # first pass: get data from all valid files
        ref = None
        for file in os.listdir(self.dir):
            fpath = self.dir+'/'+file
            try:
                assert(os.path.isfile(fpath))
                assert("median" not in file)
                assert(fpath != self.bgfile)
                hdu = fits.open(fpath)[0]
                images.append(hdu.data.astype('float'))
                header0 = hdu.header
                headers.append(header0)
            except Exception as e:
                print(red+f"{e}\nSkipping {file}"+cend)
                pass

        N_preqc = len(images)
        print(yellow+f"Pre  QC N = {N_preqc}"+cend)

        # second pass: remove corrupted images
        valid_images = []
        valid_headers = []
        nx = images[0].shape[0]
        ny = images[0].shape[1]
        for idx in np.arange(N_preqc):
            i1 = np.max([0, idx-2])
            i2 = np.min([idx+3, N_preqc])
            Nref = i2 - i1
            refimages = np.empty((nx, ny, Nref), dtype = 'float32')
            for ridx in np.arange(Nref):
                refimages[:,:,ridx] = images[ridx+i1]
            ref = np.nanmedian(refimages,axis=2)
            if self.valid(images[idx], ref, tolerance, diff_ratio) \
                and (idx+1) not in denylist:
                valid_images.append(images[idx])
                valid_headers.append(headers[idx])
            else:
                print(red+f"Skipping image {idx}"+cend)

        # free up memory before making movie
        images = valid_images
        headers = valid_headers

        Nimages = len(images)
        print(yellow+f"Post QC N = {Nimages}"+cend)

        # make movie
        fig = plt.figure()
        frames = []
        for idx in np.arange(Nimages):
            im = plt.figimage(images[idx], cmap=self.cmap, vmin = scale[0], \
                            vmax = scale[1], origin='lower', resize=True)
            if title is not None:
                plt.title(title)
            frames.append([im])
            frame = str(len(frames)).zfill(3)
            if fitsdir is not None:
                outfile = fitsdir+f"/image{frame}.fts"
                hdu_out = fits.PrimaryHDU(images[idx], headers[idx])
                hdu_out.writeto(outfile, overwrite = True)
            if framedir is not None:
                plt.savefig(framedir+f"/frame_{frame}.png")

        mov = animation.ArtistAnimation(fig, frames, interval = 50, blit = False,
              repeat = True, repeat_delay = 1000)

        print(yellow+f"Number of frames = {len(frames)}"+cend)

        mov.save(self.outfile)

