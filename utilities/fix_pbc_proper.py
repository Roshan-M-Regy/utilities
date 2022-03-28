# Roshan Mammen Regy
# roshanm.regy@tamu.edu 
# fixes pbc for each frame of a gsd trajectory
import gsd.hoomd
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
    
traj = gsd.hoomd.open(sys.argv[1])
newfile = gsd.hoomd.open(sys.argv[2],'wb')

nchains = 1000 
size = np.full(nchains,20)
print ('# of Frames: ',len(traj)) 
for i,frame in enumerate(traj):
    boxdim = frame.configuration.box[:3]
    positions = frame.particles.position
    frame.particles.positions = fix_pbc(positions,boxdim,nchains,size)
    newfile.append(frame)

newfile.close()
