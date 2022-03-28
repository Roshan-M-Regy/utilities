# Gregory L Dignon 
#!/usr/bin/python
import sys, os, numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
#import matplotlib.pyplot as plt

print ('')

traj = sys.argv[1]
top = sys.argv[2]
start = int(sys.argv[3]) # Number of frames to skip at the beginning
end = -1
#if len(sys.argv) > 4:
#    end = int(sys.argv[4])

stride = 10

u = mda.Universe(top, traj)
n_chain = len(u.atoms.segments)
N_atom = len(u.atoms.segments[0].atoms)
Rij = np.zeros((int((len(u.trajectory)-start)/stride+1), N_atom-1))
norm_factor = np.array(range(1,N_atom)[::-1])

frame = 0
for ts in u.trajectory[start:end:stride]:
    for chain in u.atoms.segments:
        dist_array = distance_array(chain.atoms.positions, chain.atoms.positions, box=u.dimensions)
        for i in range(N_atom-1):
            Rij[frame,i] = np.diagonal(dist_array, offset=i+1).mean()
    if (frame + start) % 100 == 0:
        sys.stdout.write('\rTimestep %i' % (frame + start))
        sys.stdout.flush()
    frame += 1
print ('')

## Divide into blocks and get block avg/standard deviation
n_block = 5
l_block = int(Rij.shape[0] / n_block)
block_avg = np.zeros((n_block,Rij.shape[1]))

for block in range(n_block):
    block_avg[block] = Rij[block*l_block:(block+1)*l_block].mean(0)

outp = np.zeros((N_atom-1,3))
outp[:,0] = range(1,N_atom)
outp[:,1] = block_avg.mean(0)
outp[:,2] = block_avg.std(0)/np.sqrt(n_block)

np.savetxt(sys.argv[4], outp, fmt='%-4i  %8.4f  %8.5f')
