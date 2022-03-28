# Roshan M Regy
# roshanm.regy@tamu.edu
import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis as analysis
import sys, os
from numba import jit 
import matplotlib.pyplot as plt

@jit(nopython=True)
def get_dist(x,y,z):
    # Converts distances to nm and returns the radius 
    return np.power(np.power(x/10.0,2)+np.power(y/10.0,2)+np.power(z/10.0,2),0.5)

@jit
def get_hist(subdist,nbins,boxdim,masses):
    #return np.histogram(subdist,bins=nbins,range=(0,boxdim/2.0),weights=masses)
    return np.histogram(subdist,bins=nbins,weights=masses)

@jit(nopython=True)
def get_rho(hist,bin_edges):
    for i,binno in enumerate(np.arange(len(hist))):
        #if i==0:
            #hist[i]/=4/3*np.pi*np.power(bin_edges[i+1]-bin_edges[i],3)#/1000
        #    hist[i]/=4/3*np.pi*np.power(bin_edges[i],3)#/1000
        #    continue
        hist[i] /= (4*np.pi*(np.power(bin_edges[i+1],2)*(bin_edges[i+1]-bin_edges[i])))#/1000.
        hist[i] *= 1.66053904 # Converts amu to kg/m3 
    return hist

#@jit(nopython=True)
def get_rhomat(rhomat,index,distmat,ubinedges,boxdim,masses):
    submass = masses[index-1]
    for frame in np.arange(0,distmat.shape[1],1):
        smalldist=np.zeros(len(index),dtype=np.float32)
        smalldist=distmat[index-1,frame]
        rhomat[frame,:],bin_edges = get_hist(smalldist,ubinedges,boxdim,submass) # Total mass in each bin
        rhomat[frame,:] = get_rho(rhomat[frame,:],bin_edges) # Divide mass by volume of each bin
    bins=np.zeros(len(bin_edges)-1)
    for i,val in enumerate(bin_edges[1:]):
        bins[i] = (val+bin_edges[i])/2.0
    return rhomat,bins

#@jit(nopython=True)
def get_SEM(outmat,rhomat,solconc,nblocks):
    bskip=int(rhomat.shape[0]/nblocks)
    for b in range(nblocks):
        outmat[:,b+1] = np.mean(rhomat[b*bskip:(b+1)*bskip,:],axis=0)/solconc
    outmat[:,-2] = np.mean(outmat[:,1:-2],axis=1)
    outmat[:,-1] = np.std(outmat[:,1:-2],axis=1)/np.power(nblocks,0.5)
    return outmat

def structured_binning(regions,start,stop,binno):
    binedges = [0.0]
    for i in range(regions):
        if i == 0:
            beg = start
            end = stop[0]
        else:
            beg = stop[i-1]
            end = stop[i]

        if i+1 == regions:
            nbins = binno[i]+1
            nbinedges = np.linspace(beg,end,num=nbins)
        else:
            nbins = binno[i]
            nbinedges = np.linspace(beg,end,num=nbins,endpoint=False)
        binedges = np.concatenate((binedges,nbinedges),axis=0)
    return binedges



traj = sys.argv[1]
top = sys.argv[2]
index = open(sys.argv[3], 'r').readlines()
#nbin = int(sys.argv[4])
cutoff = int(sys.argv[4])
outfolder = sys.argv[5]
noofregs = int(sys.argv[6])
start = float(sys.argv[7])
stop = np.array(sys.argv[8:8+noofregs]).astype(float)
noofbins = np.array(sys.argv[8+noofregs:8+2*noofregs]).astype(int)

blocks = 2
univ = mda.Universe(top, traj)
indexname = []
indexlist = []
for line in index:
    indexname.append(line.split()[0])
    indexlist.append(np.array(line.split()[1:]))

distmat = np.zeros((univ.atoms.positions[:,0].shape[0],len(univ.trajectory[cutoff:])))

for i,ts in enumerate(univ.trajectory[cutoff:]):
    cog = univ.atoms[indexlist[0].astype(int)-1].center_of_mass()
    #print (cog)
    x = np.array(univ.atoms.positions[:,0]-cog[0])
    y = np.array(univ.atoms.positions[:,1]-cog[1])
    z = np.array(univ.atoms.positions[:,2]-cog[2])
    distmat[:,i] = get_dist(x,y,z)

#binedges = structured_binning(2,1.0,[11.0,25.0],[50,10])
binedges = structured_binning(noofregs,start,stop,noofbins)
all_together = np.zeros((len(binedges)-1,len(indexname)))
rhomat = np.zeros((distmat.shape[1],len(binedges)-1))
for k,name in enumerate(indexname):
    print (name)
    rhomat,bins = get_rhomat(rhomat,indexlist[k].astype(int),distmat,binedges,univ.dimensions[0]/10.0,np.array(univ.atoms.masses))
    solconc = np.sum(np.array(univ.atoms.masses)[indexlist[k].astype(int)-1])/np.power(univ.dimensions[0]/10.0,3)*1.66053904
    outmat = np.zeros((len(bins),blocks+3))
    outmat[:,0] = bins
    outmat = get_SEM(outmat,rhomat,solconc,blocks)
    np.savetxt('%s/%s_mass_rho_norm_sol_conc_%s.txt'%(outfolder,traj[:9],name), outmat, header='Z  density', fmt='%10.10f')
    all_together[:,k] = outmat[:,-1]
