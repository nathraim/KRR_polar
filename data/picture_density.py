#!/usr/bin/python3
import sys
import numpy as np
import time
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.special import erf
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import LogNorm, ListedColormap
import matplotlib.pylab as pl
import matplotlib

"""
   Reads an xyz file and visualize the atomic density.
"""

name=sys.argv[1]

# Choose colormap
cmap = pl.cm.jet

# Get the colormap colors
my_cmap = cmap(np.arange(cmap.N))

# Set alpha - This way we make the grid points where the density is close to 0 transparent
#my_cmap[:,-1] = np.geomspace(5e-2, 1, cmap.N)
#my_cmap[:,-1] = np.logspace(-2, 0, cmap.N,base=10)
my_cmap[:,-1] = np.logspace(-0.2, 0, cmap.N,base=2)
#print my_cmap[:,-1]
# Create new colormap
my_cmap = ListedColormap(my_cmap)

def get_dist(point1, point2):
    a=np.array(point1)
    b=np.array(point2)
    return np.sqrt(np.sum((a-b)**2))

def density(xx,yy,zz,pos_at):
  sigma_D = 0.5 #orig
  sigma_f = 2*sigma_D
  pos_at=np.array(pos_at)
  #return np.exp(-((xx-pos_at[0])**2+(yy-pos_at[1])**2+(zz-pos_at[2])**2)/(2*sigma_D**2) - (xx**2+yy**2+zz**2)/(2*sigma_f**2)) # with filtering
  return np.exp(-((xx-pos_at[0])**2+(yy-pos_at[1])**2+(zz-pos_at[2])**2)/(2*sigma_D**2) ) # without filtering

data_geom=open(name,'r')

# Create the 3D grid
# The grid is centered at (0,0,0), which corresponds approximately to the central carbon atom

#drx,dry,drz=0.4,0.4,0.4
drx,dry,drz=0.5,0.5,0.5
xmax,ymax,zmax=6,4,0
xmin,ymin,zmin=-xmax,-ymax,-zmax
nx=int(round(2*xmax/drx+1))
ny=int(round(2*ymax/dry+1))
nz=int(round(2*zmax/drz+1))
xx, yy, zz=np.mgrid[xmin:xmax+drx:drx,ymin:ymax+dry:dry,zmin:zmax+drz:drz] #3D
SUM = xx+yy+zz
shape = SUM.shape
print ("grid size: ", SUM.shape)


#plt.ion()

trucmuche=0
while True:
  line=data_geom.readline() # Read a line that corresponds to the number of atoms
  if not line:
     break
  natoms=int(line.split()[0])
  line=data_geom.readline() # Read a second line but does nothing with it
  pos=[]
  for i in range(natoms):
     line=data_geom.readline()
     pos.append(list(map(float, line.split()[1:4]))) # x,y,z
  den=np.zeros(shape)
  for j in pos:
    den=den+density(xx,yy,zz,j)
    #den=den+density(xx,yy,0,j)
  data=[]
  for i,el1 in enumerate(den):
    for j,el2 in enumerate(el1):
      for k,el3 in enumerate(el2):
        #print xx[i,0,0],yy[0,j,0],zz[0,0,k],el3
        data.append([xx[i,0,0],yy[0,j,0],zz[0,0,k],el3])
  data=np.array(data)
        
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  
  #scat = ax.scatter(data[:,0], data[:,1], data[:,2], c=data[:,3], cmap=my_cmap, s=50, lw=0, depthshade=0)
  scat = ax.scatter(data[:,0], data[:,1], data[:,2], c=data[:,3], cmap=my_cmap, s=30, lw=0, depthshade=0)
  #plt.colorbar()
  #plt.colorbar().set_label('Atomic density', rotation=270)
  #plt.gray() 
  ax.set_xlabel('$x(\\AA$)', labelpad=8,fontsize=20,fontname='times new roman')
  ax.set_ylabel('$y(\\AA$)', fontsize=20,fontname='times new roman')
  ax.set_zlabel('$z(\\AA$)', labelpad=-5,rotation=180, fontsize=20,fontname='times new roman')
  plt.xticks(fontname = "times new roman", fontsize= 16)
  plt.yticks(fontname = "times new roman", fontsize = 16)
  ax.set_zticklabels([])
  #cb = fig.colorbar(scat,shrink=0.5,aspect=5,label='$\\rho$')
  cb = fig.colorbar(scat,shrink=0.5,aspect=5)
  axcb=cb.ax
  #text = axcb.yaxis.label
  #font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=20)
  #text.set_font_properties(font)
  axcb.text(0.25,1.05,'$\\rho$',rotation=0,fontsize=20,fontname='times new roman')
  axcb.set_yticklabels(axcb.get_yticklabels(), fontsize=16,fontname='times new roman')
  #fig.colorbar().set_label('Atomic density', rotation=270)
  plt.savefig('plot_dengrid', bbox_inches ='tight', dpi=300)
  plt.show()
    
