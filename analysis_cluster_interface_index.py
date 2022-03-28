# Roshan M Regy
# roshanm.regy@tamu.edu
# Script using freud library to calculate clusters and their interfaces with the rest of the system
import gsd.hoomd
import freud 
import numpy as np 
traj = gsd.hoomd.open('prod.gsd','rb')
# Computing clusters here
cl = freud.cluster.Cluster()
cl.compute((traj[-1].configuration.box,traj[-1].particles.position[:10000,:]),neighbors = {'r_max':1.0})
iface = freud.interface.Interface()
typeid = traj[-1].particles.typeid[10000:]
for i,val in enumerate(typeid):
    if val==0:
        typeid[i] = 1

for cluster in range(len(cl.cluster_keys)):
    # Computing interfaces with each cluster 
    iface.compute((traj[-1].configuration.box, traj[-1].particles.position[cl.cluster_keys[cluster]]),traj[-1].particles.position[10000:], neighbors={'r_max':1.5})
    print ("B particles: %s"%len(iface.query_point_ids))
    print ("A particles: %s"%len(iface.point_ids))
    indfile = open('cluster_%s_index.ndx'%(cluster+1),'w')
    AAstr = []
    ABstr = []
    indfile.write('System %s\n'%(' '.join((1+np.arange(len(cl.cluster_keys[cluster])+len(iface.query_point_ids))).astype(str))))
    indfile.write('A-A %s\n'%(' '.join((1+np.arange(len(cl.cluster_keys[cluster]))).astype(str))))
    indfile.write('X-A %s\n'%(' '.join((len(cl.cluster_keys[cluster])+1+np.arange(len(iface.query_point_ids))).astype(str))))
    X = []
    for i,val in enumerate(iface.query_point_ids):
        if traj[-1].particles.typeid[10000+val] == 1:
            X.append(len(cl.cluster_keys[cluster])+i)
    indfile.write('X %s\n'%(' '.join(np.array(X).astype(str))))
    indfile.close()
