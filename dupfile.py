"""
Since I only have 200 files for the CME model and I need 672 frames to do the 7-day noise-gate timings, this script just makes two duplicates for each file in a specified directory.
"""

import numpy as np
import os
import shutil

dir = '/home/mark.miesch/Products/image_processing/ATBD/data/timings'

Ncopy = 3

flist = os.listdir(dir)

for c in np.arange(1,Ncopy+1):
    ext = f"_{c}"
    for f in flist:
        infile = dir+'/'+f
        ff = f.split('.fts')[0] + ext + '.fts'
        outfile = dir+'/'+ff
        shutil.copy(infile,outfile)

