"""
script for making a stacked plot for the ATBD
This is for a 5-row figure made from beaf frames
"""

import imageio.v3 as iio
import numpy as np
import matplotlib.pyplot as plt

dir = '/home/mark.miesch/Products/image_processing/frames/'
outdir = '/home/mark.miesch/Products/image_processing/figs/'

fig = 6

if fig == 1:
    #file1 = dir + 'stereo_2012_09_16_ba/frame_019.png'
    file1 = dir + 'stereo_2012_09_16_ba/frame_119.png'
    file2 = dir + '2012_04_16/frame_001.png'
    outfile = outdir+'L3survey1.png'

elif fig == 2:
    file1 = dir + '2021_05_23_f6/frame_037.png'
    file2 = dir + '2021_05_noh/frame_037.png'
    outfile = outdir+'L3survey2.png'

elif fig == 3:
    file1 = dir + '2014_01_16_ba/frame_174.png'
    file2 = dir + 'CME0/frame_027.png'
    outfile = outdir+'L3survey3.png'

elif fig == 4:
    #file1 = dir + '2021_05_23_f6/frame_079.png'
    file1 = dir + '2021_05_23_f6/frame_129.png'
    #file1 = dir + '2014_01_16_ba/frame_040.png'
    file2 = dir + '2014_01_16_ba/frame_094.png'
    outfile = outdir+'L3survey4.png'

elif fig == 5:
    #file1 = dir + '2021_05_23_f6/frame_129.png'
    file1 = dir + '2014_01_16_ba/frame_051.png'
    file2 = dir + '2014_01_16_ba/frame_052.png'
    outfile = outdir+'L3survey5.png'

elif fig == 6:
    file1 = dir + 'CME0_gaussian/frame_027.png'
    file2 = dir + 'CME0_salt/frame_027.png'
    outfile = outdir+'L3noise.png'

else:
    print("enter a valid fig number")
    exit()



#------------------------------------------------------------------------------
row1 = iio.imread(file1)
row2 = iio.imread(file2)

#------------------------------------------------------------------------------

fig, ax = plt.subplots(nrows=2, ncols=1, figsize = (14,12))

plt.subplots_adjust(hspace=-0.4)

ax[0].imshow(row1)
ax[0].axis('off')

ax[1].imshow(row2)
ax[1].axis('off')

fig.tight_layout(pad=1,rect=(0.0,0.0,1.0,1.0))

xx1 = .08
yy2 = 0.92
plt.annotate("(a)", (xx1,yy2), xycoords = 'figure fraction', color='white', \
             fontsize = 'x-large', fontweight = 'semibold')

xx2 = xx1 + 0.41
plt.annotate("(b)", (xx2,yy2), xycoords = 'figure fraction', color='white', \
             fontsize = 'x-large', fontweight = 'semibold')

yy1 = 0.44
plt.annotate("(c)", (xx1,yy1), xycoords = 'figure fraction', color='white', \
             fontsize = 'x-large', fontweight = 'semibold')

plt.annotate("(d)", (xx2,yy1), xycoords = 'figure fraction', color='white', \
             fontsize = 'x-large', fontweight = 'semibold')


plt.savefig(outfile,bbox_inches='tight')


