
"""
This is a simple driver for the L3Proc class
"""
import os
import numpy as np

from astropy.io import fits
from CIMP import L3Proc as proc
from time import perf_counter

#------------------------------------------------------------------------------

dir = '/home/mark.miesch/Products/image_processing/ATBD/data'

# default parameters for all
rmin = 0.16
rmax = 1.0
clip = (0.9, 1.4)
qcrange = (0.9, 1.4)

# if Nfiles list is None, do all files in directory
Nfiles = None

# if true, wipe output directory before writing to it
# the QC filter may not function properly if you do not do this
wipe = True

fig = 2

if fig == 1:

    # L0.5 LASCO data
    # this sample of 10 includes a corrupted image, at 063005, so it 
    # should flag that if it is operating correctly
    Nfiles = 10
    endfile = dir+'/lasco_c3/L2proxy_2012_04/LASCOC3_2012_04_15_064205.fts'
    outdir = dir+'/lasco_c3/L3_2012_04'

elif fig == 2:

    # L0.5 LASCO data
    # This is the same endfile as 1 but should span 2 days
    Nfiles = 250
    endfile = dir+'/lasco_c3/L2proxyb_2012_04/LASCOC3_2012_04_16_111805.fts'
    outdir = dir+'/lasco_c3/L3b_2012_04'

elif fig == 3:

    # L0.5 LASCO data
    Nfiles = 400
    endfile = dir+'/lasco_c3/L2proxyc_2014_01/LASCOC3_2014_01_17_053005.fts'
    outdir = dir+'/lasco_c3/L3_2014_01'

elif fig == 4:

    # L1 STEREO-A data
    Nfiles = 300
    endfile = dir+'/stereo_a/L2proxy_2012_09/STEREOA_2012_09_17_022400.fts'
    outdir = dir+'/stereo_a/L3_2012_09'

elif fig == 5:

    # HAO CME model
    Nfiles = 200
    endfile = dir+'/model/CME0_pos30/L2proxy/Model0_2010_04_17_234614.fts'
    outdir = dir+'/model/CME0_pos30/L3'
    rmin = 0.1
    rmax = 1.0
    clip = (1.0, 6.0)
    qcrange = (0.0, 10.0)

elif fig == 6:

    # L0.5 LASCO data from 2021
    Nfiles = 300
    endfile = dir+'/lasco_c3/L2proxy_2021_05/LASCOC3_2021_05_23_223007.fts'
    outdir = dir+'/lasco_c3/L3_2021_05'

elif fig == 7:

    # L1 LASCO data
    Nfiles = 400
    endfile = dir+'/lasco_c3/L2proxy1_2014_01/LASCOC3_L1_2014_01_17_050449.fts'
    outdir = dir+'/lasco_c3/L3_L1_2014_01'

elif fig == 8:

    # HAO CME model with gaussian noise
    Nfiles = 200
    endfile = dir+'/model/CME0_pos30/L2proxy_gaussian/Model0_2010_04_17_234614.fts'
    outdir = dir+'/model/CME0_pos30/L3_gaussian'
    rmin = 0.1
    rmax = 1.0
    clip = (1.0, 6.0)
    qcrange = (0.0, 10.0)

elif fig == 9:

    # HAO CME model with salt noise
    Nfiles = 200
    endfile = dir+'/model/CME0_pos30/L2proxy_salt/Model0_2010_04_17_234614.fts'
    outdir = dir+'/model/CME0_pos30/L3_salt'
    rmin = 0.1
    rmax = 1.0
    clip = (1.0, 6.0)
    qcrange = (0.0, 10.0)

elif fig == 10:

    # L1 STEREO-A data - 7-day run
    Nfiles = 600
    endfile = dir+'/stereo_a/L2proxy_2012_09/STEREOA_2012_09_21_222400.fts'
    outdir = dir+'/stereo_a/L3_2012_09_7day'

else:
    print("pick a valid figure number")
    exit()

if not os.path.exists(outdir):
    os.mkdir(outdir)

#------------------------------------------------------------------------------
# clean output directory

if wipe:
    for file in os.listdir(outdir):
        fpath = outdir+'/'+file
        if os.path.isfile(fpath):
            os.remove(fpath)

#------------------------------------------------------------------------------
# get a list of files of length Nfiles, ending at endfile

indir = os.path.dirname(endfile)
dirlist = list(sorted(os.listdir(indir)))

if Nfiles == None:
    flist = dirlist
else:

    try:
        i2 = dirlist.index(os.path.basename(endfile))
    except ValueError:
        print("ERROR: file not found")
        exit()

    i1 = np.max((0,int(i2 - Nfiles + 1)))

    flist = dirlist[i1:i2+1]

#------------------------------------------------------------------------------

tstart = perf_counter()

for file in flist:
    print(80*'-')
    print(file)
    fpath = indir+'/'+file
    x = proc.l3proc(fpath, outdir)
    x.process(rmin = rmin, rmax = rmax, clip = clip, qcrange = qcrange)
    x.write()

tstop = perf_counter()

dt = tstop - tstart
dtavg = dt / len(flist)

print(f"{outdir}")
print(f"Total Time (s)   : {dt}")
print(f"Time per file (s) : {dtavg}")

#------------------------------------------------------------------------------

qc1 = []
qc2 = []
for file in os.listdir(outdir):
    hdu = fits.open(outdir+'/'+file)
    if hdu[0].header['L3QCFLAG'] == 1:
        qc1.append(file)
    if hdu[0].header['L3QCFLAG'] == 2:
        qc2.append(file)
    hdu.close()

print(80*'-'+'\n'+f"Nfiles = {Nfiles}")
print(f"QC=1   : {len(qc1)}")
for file in qc1:
    print(f"    {file}")

print(f"QC=2   : {len(qc2)}")
for file in qc2:
    print(f"    {file}")

#------------------------------------------------------------------------------
