# Roshan Mammen Regy
# roshanm.regy@tamu.edu 
# fixes pbc for each frame of a gsd trajectory
#import gsd.hoomd
import MDAnalysis as mda 
import numpy as np 
import sys
from numba import jit

@jit
def fix_pbc(positions,boxdim,nchain,size):
    start=0
    end=size[0]
    for chain in np.arange(nchain):
        for dim in np.arange(3):
            for res in np.arange(start+1,end):
                dist = positions[res-1,dim]-positions[res,dim]
                if np.absolute(dist)>=boxdim[dim]/2.0:
                    positions[res,dim] += np.sign(dist)*boxdim[dim]
        start = end
        end += size[chain+1]
    return positions
    
univ = mda.Universe(sys.argv[2],sys.argv[1])
#traj = gsd.hoomd.open(sys.argv[1])
#newfile = gsd.hoomd.open(sys.argv[2],'wb')

nchains = 70
nres = int(len(univ.atoms)/nchains)
size = np.full(nchains,nres)
print ('# of Frames: ',len(univ.trajectory)) 
w = mda.Writer(sys.argv[1][:-4]+'_pbcfix.xtc',univ.trajectory.n_atoms)
w.dt = 1000

for i,frame in enumerate(univ.trajectory):
    boxdim = univ.dimensions[:3]
    positions = univ.atoms.positions
    frame.positions = fix_pbc(positions,boxdim,nchains,size)
    w.write(univ.atoms)

w.close()

