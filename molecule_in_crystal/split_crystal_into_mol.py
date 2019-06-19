#!/usr/bin/python3
import sys

"""
This script reads an MD trajectory for paracetamol in xyz format, and breaks it into 4 trajectories for the 4 molecules contained within the unit cell 
"""

filename=sys.argv[1]

data=open(filename,'r')

prefix="mol"
Nmol=4 #Number of molecules per unit cell


for i in range(1,Nmol+1,1): # Empty previous files if they existed
  mol = open(prefix+'%i_' %i+filename, 'w')
  mol.close()

while True:
  line = data.readline()
  if not line:
    break
  Natoms = int(line.split()[0])//4 # Number of atoms per unit cell for 1 molecule
  line1 = data.readline() # Comment line
  for i in range(1,Nmol+1,1): # We loop Nmol times
    with open(prefix+'%i_' %i+filename, 'a') as f:
      f.write(str(Natoms)+'\n'+line1) # For each molecule, we print the number of atoms and 1 line of comment, to keep the same xyz file format
      for j in range(Natoms): # Prints the coordinates of current molecule
        line = data.readline()
        f.write(line)
