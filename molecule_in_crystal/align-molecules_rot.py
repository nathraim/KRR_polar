#!/usr/bin/python
import numpy as np
import sys, glob, math

"""
   Written by ???.
"""

__all__ = ['main']

def abc2h(a, b, c, alpha, beta, gamma):
   h = np.zeros((3,3) ,float)
   h[0,0] = a
   h[0,1] = b*math.cos(gamma)
   h[0,2] = c*math.cos(beta)
   h[1,1] = b*math.sin(gamma)
   h[1,2] = (b*c*math.cos(alpha) - h[0,1]*h[0,2])/h[1,1]
   h[2,2] = math.sqrt(c**2 - h[0,2]**2 - h[1,2]**2)
   return h

def invert_ut3x3(h):
   ih = np.zeros((3,3), float)
   for i in range(3):
      ih[i,i] = 1.0/h[i,i]
   ih[0,1] = -ih[0,0]*h[0,1]*ih[1,1]
   ih[1,2] = -ih[1,1]*h[1,2]*ih[2,2]
   ih[0,2] = -ih[1,2]*h[0,1]*ih[0,0] - ih[0,0]*h[0,2]*ih[2,2]
   return ih

def pbcdist(q1, q2, h, ih): 
      s = np.dot(ih,q1-q2)
      for i in range(3): s[i] -= round(s[i])
      return np.dot(h, s)


def read_frame(filedesc):
  natoms = int(filedesc.readline())
  comment = filedesc.readline()

  cell = np.zeros(6,float)
  #names=np.zeros(natoms,np.dtype('|S6')) #Nath: doesn't seem to work with Python3
  names=np.zeros(natoms,np.dtype('str'))
  q=np.zeros((natoms,3),float)
  #cell[:] = comment.split()[2:8]
  cell*=[1,1,1,0.017453293,0.017453293,0.017453293]

  for i in range(natoms):
    line=filedesc.readline().split();
    names[i]=line[0]
    q[i]=line[1:4]

  return [cell, names, q]
    
def shiftcenter(pos):
   cp=np.zeros(len(pos[0]))
   for i in range(len(pos)):
      cp+=pos[i,:]
   cp=cp/len(pos)
   spos = pos.copy()
   for v in spos: v-=cp
   #print('true centroid',cp)
   return spos

def main(prefix, reference):
   # reference structure will be shifted to zero centroid (not center of mass)
   iref=open(reference,"r")
   [ cell, names, reforg ] = read_frame(iref)
   nr = len(reforg)
   #print (nr, reforg)
   ref=shiftcenter(reforg)

# Here we read directly the position    
   ipos=open(prefix+".xyz","r")
   oposfile=open(prefix+(".align.xyz"),"w") 
# Store rotation matrix and its inverse
   orotfile=open(prefix+(".rot"),"w") 
   orotinvfile=open(prefix+("_inv.rot"),"w") 

   ifr=0
   #dt=20.670687 # 0.5 fs in atomic time. CHANGE THIS ACCORDINGLY!
   #angtobohr=1.8897261   
   #kt=np.zeros((nr,3,3),float)
   rpos=np.zeros((nr,3),float)
   rposold=np.zeros((nr,3),float)
   while True: 
     try:
       [ cell, names, posc ] = read_frame(ipos)
     except: sys.exit(0)

     nat=len(posc)
     #print nat, "\n # rotated xC"
     #if (ifr>0):
     #   ofile.write("%d\n# Velocity filtered by rotating frame. Frame: %d\n" % (nat, ifr))
     oposfile.write("%d\n# Position rotated to molecular ref. Frame: %d\n" % (nat, ifr))
     orotfile.write("# Rotation with respect to ref. Frame: %d\n" % ifr)
     orotinvfile.write("# Rotation with respect to ref. Frame: %d\n" % ifr)
     for i in range(0,nat,nr):
       #print (i, nat, nr)
       #jvel=velc[i:i+nr];
       jpos=posc[i:i+nr] 
       jnames=names[i:i+nr]
       #print 'pos', jpos
       sc, rot= kabsch(jpos, ref, names)
       #print 'sc', sc, jpos[0]-sc
       rpos[:]=jpos
       #print rpos-sc
       for v in rpos: v[:]=np.dot(rot, v-sc)
       #print 'newpos', rpos
       #Not really necessary, but move structure to zero centroid (considering all atoms), to be consistent with the ref struct
       rpos=shiftcenter(rpos)
       #Write aligned structure
       for j in range(nr): 
           oposfile.write("%6s  %15.7e  %15.7e  %15.7e\n" % (jnames[j],rpos[j,0], rpos[j,1], rpos[j,2]))
       #Write rotation matrices
       for j1 in rot:
         for j2 in j1:
             orotfile.write("%3.4e " %j2)
       orotfile.write("\n")
       invrot=np.linalg.inv(rot)
       #print 'rot', rot
       #print 'rot x invrot', np.dot(rot,rot.T)
       for j1 in invrot:
         for j2 in j1:
             orotinvfile.write("%3.4e " %j2)
       orotinvfile.write("\n")
     ref=rpos #The reference becomes the last structure
     ifr+=1

def kabsch(pos, reforg, names):

   #Calculate centroid only for the atoms considered (here, all except hydrogens)
   count=0
   cp=np.zeros(len(pos[0]))
   for i in range(len(pos)):
      if names[i]!='H':
        cp+=pos[i,:]
        count+=1
      #cp+=pos[i,:]
   #cp=cp/len(pos)
   cp=cp/count

   # now get only heavyatoms:
   indexes=[i for i,j in enumerate(names) if j!="H"]
   #indexes=[i for i,j in enumerate(names)] # We consider all atoms

   #To align only the benzene ring, uncomment the next 2 lines
   #print 'We align only the benzene ring here !'
   #indexes=indexes[4:]
   #indexes=[5, 8, 9, 11, 13, 15, 17, 18]
   spos = pos.copy()[indexes]
   ref=reforg.copy()[indexes]
   for v in spos: v-=cp

   A=np.dot(ref.T,spos)
   U, s, V = np.linalg.svd(A)      
   d = np.sign(np.linalg.det(np.dot(V.T, U.T)))
   V[2,:]*=d
   return (cp,np.dot(U, V))

if __name__ == '__main__':
   main(sys.argv[1], sys.argv[2])
