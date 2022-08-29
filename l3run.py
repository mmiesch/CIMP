
"""
This is a simple driver for the L3Proc class
"""
import os
import numpy as np

from CIMP import L3Proc as proc

#------------------------------------------------------------------------------

dir = '/home/mark.miesch/Products/image_processing/ATBD/data/'

rmin = 0.16
rmax = 1.0

fig = 1

if fig == 1:

    # L0.5 LASCO data
    infile = dir+'/lasco_c3/L2proxy_2012_04/LASCOC3_2012_04_15_063005.fts'
    outdir = dir+'/lasco_c3/L3_2012_04'

elif fig == 2:

    # L0.5 LASCO data
    infile = dir+'lasco_c3/L2proxy_2014_01/LASCOC3_2014_01_15_125405.fts'
    outdir = dir+'lasco_c3/L3_2014_01'

elif fig == 3:

    # L1 STEREO-A data
    infile = dir+'/stereo_a/L2proxy_2012_09/STEREOA_2012_09_15_183900.fts'
    outdir = dir+'/stereo_a/L3_2012_09'

elif fig == 4:

    # HAO CME model
    infile = dir+'/model/CME0_pos30/L2proxy/Model0_2010_04_16_021057.fts'
    outdir = dir+'/model/CME0_pos30/L3'
    rmin = 0.0
    rmax = np.inf


else:
    print("pick a valid figure number")
    exit()

if not os.path.exists(outdir):
    os.mkdir(outdir)

#------------------------------------------------------------------------------
#  Use the background file as a reference for the correct resolution

x = proc.l3proc_ng(infile, outdir)

#------------------------------------------------------------------------------
