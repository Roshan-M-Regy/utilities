# Roshan M Regy
# roshanm.regy@tamu.edu
import sys
lines = open(sys.argv[1]).readlines()
newfile = open(sys.argv[1][:-4]+'_CA_only.pdb','w')

for line in lines:
    if len(line.split())==0:
        continue
    if line.split()[0] not in ['ATOM','TER','END']:
        continue
    if line.split()[0]!='ATOM':
        newfile.write(line)
    else:
        if line.split()[2]=='CA':
            newfile.write(line)

newfile.close()
