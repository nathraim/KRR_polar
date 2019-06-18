#!/usr/bin/python
import sys

# This script converts a multiple line geometry to a single line one, not considering the atom types

name = sys.argv[1]

data_geom=open(name,'r')

while True:
  line=data_geom.readline()
  if not line:
    break
  line=data_geom.readline()
  for i in range(20):
    line=data_geom.readline()
    print line.split()[1],line.split()[2],line.split()[3],
  print ''
      
   
