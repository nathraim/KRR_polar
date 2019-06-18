#!/usr/bin/python3
import numpy as np
import sys

"""
   This scripts reads in a series of polarizability tensors (given as xx yy zz xy xz yz, alternating 1 line of data and 1 line of comment) and a series of (flattened) rotation matrices, and rotates the polarizabilities accordingly.
   You can choose to rotate them forward or backward.
"""

name1=sys.argv[1] # Polarizability tensors given as xx yy zz xy xz yz, alternating 1 line of data and 1 line of comment
name2=sys.argv[2] # Rotation matrices (flattened to fit on 1 line), alternating 1 line of data and 1 line of comment
name3=sys.argv[3] # Option: Forward 'f' or backward 'b' transformation

infile1 = open(name1,'r')
infile2 = open(name2,'r')

polarout = np.zeros((3,3))
counter = 0

while(True):
  counter+=1
  line1 = infile1.readline() # Pass 1 line of comment
  if (not line1):
    break
  line2 = infile2.readline() # Pass 1 line of comment
  if (not line2):
    break
  line1 = infile1.readline()
  line2 = infile2.readline()
  polarin = list(map(float,line1.split()[0:6]))
  rot = list(map(float,line2.split()[0:]))
  # Recast rotation matrix as a....matrix
  rot = np.reshape(rot,(3,3))
  # Recast the polarizability into tensorial form
  polarout[0,0] = polarin[0]
  polarout[0,1] = polarin[3]
  polarout[0,2] = polarin[4]
  polarout[1,0] = polarin[3]
  polarout[1,1] = polarin[1]
  polarout[1,2] = polarin[5]
  polarout[2,0] = polarin[4]
  polarout[2,1] = polarin[5]
  polarout[2,2] = polarin[2]

  # Now apply forward or backward transformation
  # Note that the transpose of a rotation matrix is equal to its inverse
  if  (name3 == 'f'):
    polarout = np.dot(np.dot(rot,polarout),rot.T)  # R alpha R^{-1}
  elif (name3 == 'b'):
    polarout = np.dot(np.dot(rot.T,polarout),rot)  # R^{-1} alpha R 
  print ("#Step", counter)
  print (polarout[0,0],polarout[1,1],polarout[2,2],polarout[0,1],polarout[0,2],polarout[1,2])
