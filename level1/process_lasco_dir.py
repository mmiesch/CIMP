
import glob
import numpy as np
import os
import subprocess

from astropy.io import fits
from astropy.time import Time

"""
This function is for processing STEREO-A total brightness (tB) files in "batch" mode from L0 to L1 using the SolarSoft `secchi_prep` routine.  It processes all files in a specified target directory and writes the output as fits files to a specified output directory.

"""
#------------------------------------------------------------------------------
red = '\033[91m'
yellow = '\033[93m'
cend = '\033[0m'

#------------------------------------------------------------------------------
# define platform-specific parameters

# this is the input directory
#targetdir = '/home/mark.miesch/data/lasco_monthly/c3/2014_01/'
targetdir = '/home/mark.miesch/data/lasco_monthly/c3/2021_05/'

#output directory (+ month/day)
#outdir = '/home/mark.miesch/data/lasco_monthly/c3/L1/2014_01/'
outdir = '/home/mark.miesch/data/lasco_monthly/c3/L1/2021_05/'

# location of sswidl executable
sswidl = "/usr/local/ssw/gen/setup/ssw_idl"

#------------------------------------------------------------------------------
# check input specification

if targetdir[-1] != '/':
    print(red+"Error: target directory must end in '/'"+cend)
    exit(1)

if outdir[-1] != '/':
    print(red+"Error: output directory must end in '/'"+cend)
    exit(1)

if not os.path.exists(outdir):
    os.mkdir(outdir)

#------------------------------------------------------------------------------
# loop over month directory

for subdir in os.listdir(targetdir):
    ddir = targetdir+'/'+subdir
    if os.path.isdir(ddir):
        day = subdir
        print(80*'-'+f"\nday {day}")

        savepath = outdir + '/' + day
        if not os.path.exists(savepath):
            os.makedirs(savepath)

        for file in os.listdir(ddir):
            fpath = ddir+'/'+file
            try:
                assert("median" not in file)
                assert("fts" in file)
                print(file)
                hdu = fits.open(fpath)
                nx = hdu[0].header['NAXIS1']
                ny = hdu[0].header['NAXIS2']
                assert(nx == 1024)
                assert(ny == 1024)
            except:
                print(red+f"Skipping file {file}"+cend)
                continue

            idlcommand = f"reduce_level_1,'{fpath}',savedir='{savepath}'"
            subprocess.run([sswidl,"-e",idlcommand], env=os.environ)

