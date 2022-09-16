
"""
This loops over a directory of L3 files and tells which have a nonzero QC flag
"""
import os
import numpy as np

from astropy.io import fits
from CIMP import L3Proc as proc
from time import perf_counter

#------------------------------------------------------------------------------


dir = '/home/mark.miesch/Products/image_processing/ATBD/data/lasco_c3/L3_2014_01'

#------------------------------------------------------------------------------

qc1 = []
qc2 = []
for file in os.listdir(dir):
    hdu = fits.open(dir+'/'+file)
    if hdu[0].header['L3QCFLAG'] == 1:
        qc1.append(file)
    if hdu[0].header['L3QCFLAG'] == 2:
        qc2.append(file)
    hdu.close()

print(80*'-'+'\n')
print(f"QC=1   : {len(qc1)}")
for file in qc1:
    print(f"    {file}")

print(f"QC=2   : {len(qc2)}")
for file in qc2:
    print(f"    {file}")

#------------------------------------------------------------------------------
