#!/usr/bin/python3
import sys
import numpy as np

"""
   This script converts regular text data (here assuming 1 line of comment and 1 line of data) into numpy format.
   Numpy format is generally less heavy, and the reading time can be 100-1000 times faster.
"""

name0 = sys.argv[1]
file0 = open(name0,'r')

array = []
while (True):
  line = file0.readline() # Discard comment line
  if not line:
    break
  line = file0.readline() # Line containing data
  array.append(list(map(float, line.split()[0:])))

file0.close()

array = np.array(array)
# Save to file
np.save(name0+'.npy',array)

