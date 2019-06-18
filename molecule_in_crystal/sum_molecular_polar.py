#!/usr/bin/python3
import sys
import numpy as np

"""
    This script simply reads the polarizabilities of different molecules and sum them up 
"""

Nmol = 4 # Number of molecules per unit cell
#Nconf = 15091 # Number of configurations per file
Nconf = 480 # Number of configurations per file
poltot = np.zeros((Nconf,6)) # Sum of molecular polarizabilities

for i in range(1,Nmol+1,1): # We loop Nmol times
  with open('polar_harmo_mol%i.align' %i+'_rotasincrystal.dat', 'r') as f_:
    data = f_.readlines()
    for i in range(1,len(data)+1,2):
      # Add current molecular polarizability to the total
      poltot[(i-1)/2]+=np.array(map(float,data[i].split()[0:6]))

# Print out the results
for j,el1 in enumerate(poltot):
  print ('#Step',j+1)
  for el2 in el1:
    print (el2),
  print ('\n'),

