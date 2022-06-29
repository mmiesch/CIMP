# CIMP

CIMP = Coronagraph Image Processing
CCOR = Compact Coronagraph, developed by the US Naval Research Laboratory (NRL)
SWFO-L1 = Space Weather Follow On - L1 Lagrange point

This code repository contains experimental and prototype python software developed for use in assessing and producing Level 3 data products for the CCCOR-1 instrument on NOAA's GOES-U spacecraft and the CCOR-2 instrument on NOAA's SWFO-L1 spacecraft.

The CIMP software package is located in the `CIMP` subdirectory.  The root directory contains a variety of driver routines for plotting and analysing data using the CIMP software package.  As such, the files in the root directory are less organized and more customizable.


---
## Example usage

To load the module, make sure this directory is in your path and then do one of the following:

```python
import CIMP
x = CIMP.Event.event()
```

```python
from CIMP import Event as ev
x = ev.event()
```

See the driver routines in the root directory for other examples.