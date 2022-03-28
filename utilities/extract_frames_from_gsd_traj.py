# Roshan M Regy
# roshanm.regy@tamu.edu
import gsd.hoomd
import sys

traj = gsd.hoomd.open(sys.argv[1])

start = int(sys.argv[2])

if start>len(traj):
    print ("Start is greater than number of frames %s."%len(traj))
    exit()

stop = int(sys.argv[3])
if stop>len(traj):
    print ("Stop is now equal to %s"%len(traj))
    stop=len(traj)

skip = int(sys.argv[4])
nfilename = sys.argv[5]
newfile = gsd.hoomd.open(nfilename,'wb')

for frame in traj[start:stop:skip]:
        newfile.append(frame)
newfile.close()
