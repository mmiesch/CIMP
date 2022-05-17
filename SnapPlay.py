"""
This is similar to EnhancePlay.py but it is a simplification, focusing on the development of the CIMP Snapshot class and the enhancement methods it calls through the CIMP Enhance class.  So, use EnhancePlay.py to assess new possible enhancement algorithms.  But, use this to test the algorithms that are already implemented in the Enhance class.
"""

import numpy as np
import matplotlib.pyplot as plt

from CIMP import Snapshot as snap

pcase = 1

if pcase == 1:
    testcase = 1
    scale = None
else:
    print("specify a valid test case")
    exit()

#------------------------------------------------------------------------------

x = snap.snapshot.testcase(testcase)

