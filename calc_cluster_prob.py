import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-traj", type=str, required=True,help="gsd trajectory")
parser.add_argument("-nchain",type=int, required=True,help="# of chains")
parser.add_argument("-start",type=int, required=True,help="start frame")
parser.add_argument("-skip",type=int, required=True,help="skip frames")
parser.add_argument("-out",type=str,required=True,help="Name of output file")
parser.add_argument("-blocks",type=int,required=True,help="Number of blocks for SEM")
args = parser.parse_args()
import gsd.hoomd
import freud
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_SEM(data,blocks):
    blocksize = int(data.shape[0]/blocks)
    if len(data.shape)>1:
        blockavg = np.zeros((blocks,data.shape[1]))
    else:
        blockavg = np.zeros(blocks)
    for b in range(blocks):
        if len(data.shape)>1:
            blockavg[b,:] = np.mean(data[blocksize*b:blocksize*(b+1),:],axis=0)
        else:
            blockavg[b] = np.mean(data[blocksize*b:blocksize*(b+1)])
    
    return np.mean(blockavg,axis=0),np.std(blockavg,axis=0)/np.sqrt(blocks)



traj = gsd.hoomd.open(args.traj)
framelist = range(args.start,len(traj),args.skip)
cl = freud.cluster.Cluster()
nchain = args.nchain
nres = int(traj[0].particles.N/nchain)
clustsizes = np.zeros((len(framelist),nchain))

for i,f in enumerate(framelist):
    frame =traj[f]
    cl.compute((frame.configuration.box,frame.particles.position),neighbors={'r_max':7.5})
    for c,clust in enumerate(cl.cluster_keys):
        clustsizes[i,int(len(clust)/nres)-1] +=1
    for ch,val in enumerate(clustsizes[i,:]):
        clustsizes[i,ch] = clustsizes[i,ch]/nchain*(ch+1)

clustsizedist, err = get_SEM(clustsizes, args.blocks)
outputmat = np.zeros((nchain,3))
outputmat[:,0] = range(1,nchain+1)
outputmat[:,1] = clustsizedist
outputmat[:,2] = err
df = pd.DataFrame(outputmat,columns=['Nc','P(Nc)','P(Nc) SEM'])
df.to_csv(args.out+'.csv')
fig,ax = plt.subplots(1,1,dpi=300,figsize=[3,3])
ax.errorbar(range(nchain),clustsizedist,yerr=err,marker='o',lw=0.5,markersize=2)
ax.set_xticks(range(nchain)[::20])
ax.set_xticklabels(range(1,nchain+1)[::20])
ax.tick_params(labelsize=5)
ax.set_xlabel(r'Cluster size $(N_{c})$',fontsize=7)
ax.set_ylabel(r'P$(N_{c})$',fontsize=7)
ax.set_xscale('log')
ax.set_ylim(-0.1,1)
plt.tight_layout()
plt.savefig(args.out+'.png',dpi=300)
plt.show()
