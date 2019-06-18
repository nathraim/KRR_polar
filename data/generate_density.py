#!/usr/bin/python
import sys
import numpy as np

"""
   Reads a trajectory file (xyz format) and computes the atomic density on a grid.
   The important parameters to change manually are x-y-zmax-min (extension of your system), and drx-y-z (grid spacing).
   Note: the scaling could be largely improved.
"""

name = sys.argv[1]

def get_dist(point1, point2):
    a = np.array(point1)
    b = np.array(point2)
    return np.sqrt(np.sum((a-b)**2))

def density(xx,yy,zz,pos_at):
  sigma_D = 0.5 #orig, in Angstroms
  #sigma_D = 0.7
  sigma_f = 2*sigma_D
  pos_at = np.array(pos_at)
  return np.exp(-((xx-pos_at[0])**2+(yy-pos_at[1])**2+(zz-pos_at[2])**2)/(2*sigma_D**2) ) # without filtering

data_geom = open(name,'r')

# Create the 3D grid
# The grid is centered at (0,0,0), which corresponds approximately to the central carbon atom

#drx,dry,drz=0.5,0.5,0.5 # Usual
#xmax,ymax,zmax=6,4,2.5
drx,dry,drz=0.5,0.5,0.5
xmax,ymax,zmax=1,2,2
xmin,ymin,zmin=-xmax,-ymax,-zmax
nx=int(round(2*xmax/drx+1))
ny=int(round(2*ymax/dry+1))
nz=int(round(2*zmax/drz+1))
xx, yy, zz=np.mgrid[xmin:xmax+drx:drx,ymin:ymax+dry:dry,zmin:zmax+drz:drz]
SUM=xx+yy+zz
#print 'SUM', SUM
shape=SUM.shape
print shape,
dim = shape[0]*shape[1]*shape[2]
print dim

densave = []
counter = 0
while True:
  counter+=1
  line = data_geom.readline() # Reads a line that corresponds to the number of atoms
  if not line:
     break
  natoms = int(line.split()[0])
  line = data_geom.readline() # Reads a second line but does nothing with it
  pos = []
  for i in range(natoms):
     line = data_geom.readline()
     pos.append(map(float, line.split()[1:4])) # x,y,z
  den = np.zeros(shape)
  for j in pos:
    den = den+density(xx,yy,zz,j)
   
  denflat = den.flatten() # Recasts den as a 1D array (easier to iterate over the elements)
  #denflat[denflat<1e-8] = 0 # Replace small elements by 0 to avoid writing big numbers in the output for nothing, it will thus reduce the size of it

  print counter # Just to follow the process

  #fichier.write('#Step '+str(counter)+'\n')
  #for elem in denflat:
  #  fichier.write( ("{:.4e}".format(elem)) + ' ')
  #fichier.write('\n')

  #Divide by the number of grid points to kind of normalize (this way the hyperparameters won't change magnitude if we increase or decrease the number of grid points
  denflat = denflat/dim
  densave.append(denflat)

#Save as numpy array
np.save("density_" + name + ".npy", densave)
