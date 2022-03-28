# Roshan Mammen Regy
# roshanm.regy@tamu.edu 
# fixes pbc for each frame of a gsd trajectory
import gsd.hoomd
import numpy as np 
import sys
from numba import jit

@jit
def get_res_to_fix(images):
    corrections = images-images[0]
    return corrections
@jit
def fix_pbc(positions,corrections,boxdim):
    for res in np.arange(positions.shape[0]):
        positions[res,:] += corrections[res]*boxdim[:3]
    return positions
@jit
def calling_func(nchains,images,positions,boxdim):
    for chain in np.arange(nchains):
        chainindex = np.arange(size*chain,size*(chain+1))
        cimages = images[chainindex]
        cpositions = positions[chainindex]
        corrections = get_res_to_fix(cimages)
        positions[chainindex] = fix_pbc(cpositions,corrections,boxdim)
    return positions
    
traj = gsd.hoomd.open(sys.argv[1])
newfile = gsd.hoomd.open(sys.argv[2],'wb')
# Only for rigid body only simulations. For others use bond information
bodytypes, counts = np.unique(traj[0].particles.body, return_counts=True)
nchains = len(bodytypes)
size = counts[0]
print ('# of Frames: ',len(traj)) 
for i,frame in enumerate(traj):
    print (i+1)
    boxdim = frame.configuration.box
    images = frame.particles.image
    positions = frame.particles.position
    frame.particles.positions = calling_func(nchains,images,positions,boxdim)
    newfile.append(frame)

newfile.close()
