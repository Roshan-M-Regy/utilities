# Author: Roshan Mammen Regy
# EmailID: roshanm.regy@tamu.edu
import numpy as np 
import gsd.hoomd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import sys
import string
from numba import jit

# Function to check for contact formation
# Distance based criteria, <=1.5*sigma 
# Currently sigma = 5 so the value is hardcoded 
@jit
def get_dist(arr1, arr2):
    dist = 0
    for i in range(arr1.shape[0]):
        for j in range(arr2.shape[0]):
            if (np.sum(np.power(arr1[i,:]-arr2[j,:],2))) <= 56.25:
                dist +=1
    return dist

# Calculates chain level contact matrix for passed positions
@jit
def get_nu_cmat(nchain,chainlen,positions,headpos,tailpos,boxsize,cmat):
    for c1 in (range(nchain-1)):
        for c2 in range(c1+1,nchain):
            shifts = np.array([[0,0,0],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
            htdist = 0
            thdist = 0
            for img in range(shifts.shape[0]):
                htdist += get_dist(positions[chainlen*c1+headpos],positions[chainlen*c2+tailpos]+shifts[img,:]*boxsize)
                thdist += get_dist(positions[chainlen*c1+tailpos],positions[chainlen*c2+headpos]+shifts[img,:]*boxsize)
            cmat[c1,c2]  = htdist  
            cmat[c2,c1]  = -1*thdist
    return cmat

# Builds a tree of oligomers from the contact matrix 
def get_tree(c1,cmat,j,temp_cluster):
    nxht = np.where(cmat[c1,:] !=0)[0]
    nxth = np.where(cmat[:,c1] !=0)[0]
    temp = []
    temp.extend(list(nxht))
    temp.extend(list(nxth))
    for item in list(nxht):
        if item not in temp_cluster:
            temp_cluster.append(item)
    for item in list(nxth):
        if item not in temp_cluster:
            temp_cluster.append(item)

    for i in temp_cluster:
        if i not in tree:
            tree.append(i)
            get_tree(i,cmat,j, temp_cluster)
            if temp == 0:
                continue
    cluster_list.append(temp_cluster)
    return temp_cluster

# Call this function at the end of main loop to plot chain lvl contact matrix for a frame
def plot_cmat(cmat):
    fig,ax = plt.subplots(1,1,dpi=300,figsize=[3,3])
    ticks = list(string.ascii_uppercase)+['a','b','c','d']
    ax.set_xticks(np.arange(-.5, 30, 1), minor=True)
    ax.set_yticks(np.arange(-.5, 30, 1), minor=True)
    ax.tick_params(which='minor',width=0)
    norm = colors.TwoSlopeNorm(vmin=-14, vcenter=0, vmax=14)
    ax.grid(which='minor', color='black', linestyle='-', linewidth=0.5)
    ax.set_xticks(range(30))
    ax.set_xticklabels(ticks,fontsize=4)
    ax.set_yticks(range(30))
    ax.set_yticklabels(ticks,fontsize=4)
    im = ax.imshow(cmat,origin='lower',cmap=cm.seismic,norm=norm)
    cb = fig.colorbar(im,shrink=0.5)
    cb.ax.tick_params(labelsize=5)
    ax.set_title('Frame %s'%frame,fontsize=7)
    plt.tight_layout()
    plt.show()

# Load trajectory
traj = gsd.hoomd.open(sys.argv[1])
# Scan parameters
scanparam="15 10 5 4.5 4 3.5 3 2.5 2 1.5 1 0.99 0.98 0.97 0.96 0.90 0.85 0.80 0.75 0.70 0.65 0.60 0.55 0.50 0.45".split()
# Chain properties
chainlen = 178
nchain = int(traj[0].particles.N/chainlen)
headpos = [] # Positions of the head interface
tailpos = [] # Positions of the tail interface 
# Replace above two lists with all positions if no specific interfaces are present
for i,typ in enumerate(traj[0].particles.typeid[:chainlen]):
    if typ==3:
        headpos.append(i)
    if typ==4:
        tailpos.append(i)
headpos = np.array(headpos)
tailpos = np.array(tailpos)
# Start and skip frames for calculation
start = 49
skip = 30
# Output matrix containing size and number of clusters
outputmat = np.zeros((len(scanparam),5)) # Biggest oli, Big oli std., # of oli, # of oli std. 
# Main loop 
for i,lam in enumerate(scanparam):
    print ('##############################')
    print ('ScanParam ',scanparam[i])
    startf = start+i*skip
    endf = startf+10
    print (startf,endf)
    foldmat = np.zeros((10,2)) # Biggest oli, # of oli
    for f,frame in enumerate(range(startf,endf)):
        print ('FRAME ',frame)
        positions = traj[frame].particles.position
        boxsize = traj[frame].configuration.box[:3]
        cmat = np.zeros((nchain,nchain))
        cmat = get_nu_cmat(nchain,chainlen,positions,headpos,tailpos,boxsize,cmat)
        clustermat = np.zeros(nchain)
        clustermatid = 0
        for chainid in range(nchain):
            temp_cluster = []
            global tree 
            tree = []
            global cluster_list 
            cluster_list = []
            if clustermat[chainid] ==0:
                cl = get_tree(chainid,cmat,i, temp_cluster) 
                if len(cl)==0:
                    continue
                clustermatid +=1
                clustermat[cl] = clustermatid
        types,counts = np.unique(clustermat,return_counts=True)
        mx = 0
        for t,typ in enumerate(types):
            if typ!=0:
                if counts[t]>mx:
                    mx=counts[t]
        foldmat[f,0] = mx
        foldmat[f,1] = clustermatid
    
    print (np.mean(foldmat[:,0]))
    outputmat[i,1] = np.mean(foldmat[:,0])
    outputmat[i,2] = np.std(foldmat[:,0])
    outputmat[i,3] = np.mean(foldmat[:,1])
    outputmat[i,4] = np.std(foldmat[:,1]) 

outputmat[:,0] = np.array(scanparam).astype('float')
np.savetxt('cluster_size_number.txt',outputmat,header='# Scanparam Largest_cluster Deviation #_of_clusters Deviation') 
# Plots biggest cluster size wrt scan parameter 
fig,ax = plt.subplots(1,1,dpi=300,figsize=[3,3])
ax.errorbar(-1*outputmat[:,0],outputmat[:,1],yerr=outputmat[:,2],marker='o',markersize=2,lw=0.5)
ax.set_xticks(-1*outputmat[:,0])
ax.set_xticklabels(outputmat[:,0])
plt.tight_layout()
plt.show()
#############################
