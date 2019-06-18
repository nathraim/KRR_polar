#!/usr/bin/python
import numpy as np
import sys


"""
   This scripts reads in a series of geometries (xyz format) and a series of (flattened) (inverse) rotation matrices, and rotates the coordinates accordingly.
"""

name1=sys.argv[1] # geometries
name2=sys.argv[2] # Inverse rotation matrices

infile1=open(name1,'r')
infile2=open(name2,'r')

#rot=np.array([[1,0,0],[0,-1,0],[0,0,-1]])
#rot=np.array([[-4.9148e-02, 5.7615e-01, 8.1586e-01], [-9.8151e-01, 1.2342e-01, -1.4628e-01], [-1.8497e-01, -8.0797e-01, 5.5944e-01]])
#print rot

counter=0
while(True):
  counter+=1
  line1 = infile1.readline()
  if (not line1):
    break
  natoms = int(line1.split()[0])
  print line1,
  line1 = infile1.readline()
  line2 = infile2.readline()
  print line1,
  line2 = infile2.readline()
  rot = map(float,line2.split()[0:])
  rot = np.reshape(rot,(3,3))
  for i in xrange(natoms):
    line1 = infile1.readline()
    attyp = line1.split()[0]
    geom = np.array(map(float,line1.split()[1:4]))
    newgeom = np.dot(rot.T,geom)
    print attyp, newgeom[0], newgeom[1], newgeom[2]
