"""
A tool to plot out an event
"""

from CIMP import Event as ev

x = ev.event.testcase(1)

#x = ev.event.fromtime()

print(80*'-')
print(x)
print(80*'-')
print(repr(x))
print(80*'-')


# ===================
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u

fig = plt.figure(figsize=[18,8])

for i in np.arange(0,5):
    amap = x.map(i)
    ax = fig.add_subplot(2,3,i+1,projection=amap)
    amap.plot(clip_interval=[10,90]*u.percent)

plt.show()