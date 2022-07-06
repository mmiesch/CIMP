"""
Figure illustrating morphological image processing
based on the scikit-image examples here:

https://scikit-image.org/docs/dev/auto_examples/applications/plot_morphology.html
"""

import matplotlib.pyplot as plt
from skimage import data
from skimage.morphology import disk, erosion, dilation
from skimage.util import img_as_ubyte

orig_phantom = img_as_ubyte(data.shepp_logan_phantom())

#------------------------------------------------------------------------------
# erosion and dilatation

footprint = disk(6)

eroded = erosion(orig_phantom, footprint)

dilated = dilation(orig_phantom, footprint)


#------------------------------------------------------------------------------


fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(12,4), sharex = True, \
                                    sharey = True)


ax1.imshow(orig_phantom, cmap=plt.cm.gray)
ax1.set_title('Original')
ax1.axis('off')

ax2.imshow(eroded, cmap=plt.cm.gray)
ax2.set_title('Erosion')
ax2.axis('off')

ax3.imshow(dilated, cmap=plt.cm.gray)
ax3.set_title('Dilation')
ax3.axis('off')

fig.tight_layout(pad=1,rect=(0.01,0.01,.99,.98))

yy = 0.84
dx = 0.02

print(f"MSM {ax1.get_position().x0}")

xx = ax1.get_position().x0 + dx
plt.annotate("(a)", (xx,yy), xycoords = 'figure fraction', color='white', \
             fontsize = 'x-large', fontweight = 'semibold')

xx = ax2.get_position().x0 + dx
plt.annotate("(b)", (xx,yy), xycoords = 'figure fraction', color='white', \
             fontsize = 'x-large', fontweight = 'semibold')

xx = ax3.get_position().x0 + dx
plt.annotate("(c)", (xx,yy), xycoords = 'figure fraction', color='white', \
             fontsize = 'x-large', fontweight = 'semibold')

plt.show()