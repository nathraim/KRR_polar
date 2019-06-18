#!/usr/bin/python3
import sys
"""
   Just add a line of comment on each line
"""

name = sys.argv[1]

data_geom = open(name,'r')
step = 1

while True:
  line = data_geom.readline()
  if not line:
    break
  print ("#Step", step)
  print (line.split('\n')[0])
  step += 1
