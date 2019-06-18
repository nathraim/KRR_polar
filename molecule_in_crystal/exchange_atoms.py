#!/usr/bin/python3
import sys
import numpy as np

"""
   This scripts reads an MD trajectory for the paracetamol molecule, and swaps 4 carbons within the benzene ring and their 4 hydrogens.
"""

filename = sys.argv[1]

# Read data
with open(filename, 'r') as file0:
  data = file0.readlines()

natoms = int(data[0].split()[0]) # Number of atoms
Nsteps = len(data)//(natoms+2) # Number of steps in the file

i = 1 # Index (will correspond to the numbering of the atoms)

for N in range(0,Nsteps,1): 
  data_old = data[i+10:i+14]
  data[i+10:i+14] = data[i+14:i+18]
  data[i+14:i+18] = data_old

  i+=natoms+2 # Go to the next structure (each step contains natoms+2 lines)

# Write everything back
with open('exchanged_'+filename, 'w') as file0:
  file0.writelines( data )
