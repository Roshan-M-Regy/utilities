#!/usr/bin/env python

import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis as ana
import MDAnalysis.analysis.distances 
#import time
import sys

def CalcCutoffMatrix(seq):
    cutoff = np.zeros((len(seq),len(seq)))
    for i in range(len(seq)):
        for j in range(i,len(seq)):
            cutoff[i,j] = (siglist[AAlist.index(seq[i])] + siglist[AAlist.index(seq[j])])/2.
            cutoff[j,i] = (siglist[AAlist.index(seq[i])] + siglist[AAlist.index(seq[j])])/2.
    return cutoff

# input xtc file
traj = sys.argv[1]
top = sys.argv[2]
fout = sys.argv[3]
#molname = sys.argv[4]
seqfile = sys.argv[4]
start = int(sys.argv[5])
stride = int(sys.argv[6])
stop = int(sys.argv[7])

# 
global siglist, AAlist
siglist = np.loadtxt('/raid/data/idp/move/nij218/VDW_params.dat', usecols=1)
AAlist = list('ACDEFGHIKLMNPQRSTVWY')

# input gro or pdb file
univ = mda.Universe(top,traj)

natoms = len(univ.atoms)
nchain = 2
nres = int(natoms/nchain)
#cutoff = np.loadtxt('/raid/data/idp/gld215/LAMMPS/Codes/cont_dist/%s.dat'%(molname))*(2**(1./6.))
#cutoff = CalcCutoffMatrix(open(seqfile, 'r').read().strip()) * (2**(1./6.))
cutoff = CalcCutoffMatrix(open(seqfile, 'r').read().strip()) * 1.5
style = "serial"


frame = start
# trajectory length, number of frames
nstep = int((len(univ.trajectory)-start)/stride + 1)
print (len(univ.trajectory))
# output matrix
countall = np.zeros((nstep, nres, nres),dtype=np.uint8)
for ts in univ.trajectory:
    print ('start')
    for i in range(nchain):
        for j in range(i+1,nchain):
            moli = univ.atoms.positions[i*nres:(i+1)*nres,:]
            molj = univ.atoms.positions[j*nres:(j+1)*nres,:]
            d = ana.distances.distance_array(moli, molj, box=univ.dimensions, backend=style)
            print (d)
            countall[int((frame-start)/stride)] += (d <= cutoff)
    frame += 1*stride

np.save(fout,countall)
