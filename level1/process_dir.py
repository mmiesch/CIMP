
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

# this is the directory to monitor for new files
targetdir = '/home/mark.miesch/sunpy/data/secchi_cor2/2012/'

#output directory (+ month/day)
outdir = '/home/mark.miesch/sunpy/data/secchi_cor2/L1/2012/'


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
# first get the file list of the target directory

files = list(filter(os.path.isfile, glob.glob(targetdir + "*.fts")))
files.sort(key=lambda x: os.path.getmtime(x))

#------------------------------------------------------------------------------
# now loop over files, looking for a complete set

# remaining files to be processed
filelist = files.copy()

for file in files:

    fname = file.split('/')[-1]

    hdu = fits.open(file)[0]
    time = Time(hdu.header['DATE-OBS'])

    month = time.strftime('%m')
    day = time.strftime('%d')

    nx = hdu.header['NAXIS1']
    ny = hdu.header['NAXIS2']

    savepath = outdir + month + '/' + day

    if (nx == 2048) and (ny == 2048):

        print(yellow+f"{fname} {month} {day} {nx} {ny}"+cend)

        if not os.path.exists(savepath):
            os.makedirs(savepath)


    idlcommand = f"secchi_prep,'{file}',savepath='{savepath}',/write_fts"
    subprocess.run([sswidl,"-e",idlcommand], env=os.environ)

