"""
Use this to 
"""

import numpy as np
from CIMP import Event as ev
from sunpy.net import attrs as a
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map

#---------------------------------------------------

tcase = 1

if tcase == 1:
    testcase = 1
    plotframe = 9
    scale = (0.0, 1000.0)

elif tcase == 2:
    testcase = 8
    plotframe = 3
    scale = (0.0, 1000.0)

else:
    print("specify a valid tcase")
    exit()    

#---------------------------------------------------
# Define base image

x = ev.event.testcase(testcase)

a = x[plotframe]
amap = x.map(plotframe)

header = x.header[plotframe]

print(80*'-')
print(x)
print(80*'-')
print(repr(x))
print(80*'-')

print(f"Data range: {amap.min()} {amap.max()}")

#---------------------------------------------------
# Choose your battle

comp = (0,2)

images = []
scales = []
titles = []

if comp.count(0) > 0:
    titles.append("Base image")
    images.append(a)
    scales.append(scale)

if comp.count(1) > 0:
   titles.append("NRGF")
   x.nrgf()
   b = x[plotframe]
   images.append(b)
   scales.append((0.0,4.0))

if comp.count(2) > 0:
    # current enhance() method in Event class
    titles.append("Event.enhance()")
    x.enhance(clip = scale)
    b = x[plotframe]
    images.append(b)
    scales.append(None)

# ===================

fig = plt.figure(figsize=[16,8])

map1 = sunpy.map.Map(images[0],header)
ax = fig.add_subplot(1,2,1,projection=map1)
if scales[0] is None:
    map1.plot(title=titles[0])
else:
   print(f"image 1 scale: {scales[0][0]} to {scales[0][1]}")
   map1.plot(vmin=scales[0][0], vmax=scales[0][1], title=titles[0])

map2 = sunpy.map.Map(images[1],header)
ax = fig.add_subplot(1,2,2,projection=map2)
if scales[1] is None:
    map2.plot(title=titles[1])
else:
   print(f"image 2 scale: {scales[1][0]} to {scales[1][1]}")
   map2.plot(vmin=scales[1][0], vmax=scales[1][1], title=titles[1])

# ===================
# save to a file
# ===================
dir = '/home/mark.miesch/Products/image_processing/CIMP/images/Enhance_comp/'

file = dir + f"Enhance_t{testcase}_{comp[0]}_vs_{comp[1]}.png"

plt.savefig(file)

plt.show()