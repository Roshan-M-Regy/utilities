import gsd.hoomd
import sys 
import os 
import string
import numpy as np
amino_1to3 = {
        'A':'ALA',
        'C':'CYS',
        'D':'ASP',
        'E':'GLU',
        'F':'PHE',
        'G':'GLY',
        'H':'HIS',
        'I':'ILE',
        'K':'LYS',
        'L':'LEU',
        'M':'MET',
        'N':'ASN',
        'P':'PRO', 
        'Q':'GLN',
        'R':'ARG',
        'S':'SER',
        'T':'THR',
        'V':'VAL',
        'W':'TRP',
        'Y':'TYR',
        }

traj = gsd.hoomd.open(sys.argv[1],'rb')
seq = list(open(sys.argv[2]).read().strip())
nres = len(seq)
nchain =   int(traj[0].particles.N/nres) #int(sys.argv[3])
letters = list(string.ascii_uppercase)
root = sys.argv[1][:-4]+'2pdb'
cnt = 0
outdir= root
if os.path.isdir(outdir) == False:
    os.mkdir(outdir)
for f,frame in enumerate(traj):
    outfile = open(outdir+'/frame_%s.pdb'%f,'w')
    for chain in range(nchain):
        for part in range(nres):
            spacing = [7,4,5,2,4,12,8,8]
            lineelements = ['ATOM',chain*nres+part+1,'CA', amino_1to3[seq[part]], letters[int(frame.particles.types[frame.particles.typeid[chain*nres+part]])],part+1,np.around(frame.particles.position[chain*nres+part,0],decimals=3),np.around(frame.particles.position[chain*nres+part,1],decimals=3),np.around(frame.particles.position[chain*nres+part,2],decimals=3)]
            outstring = ''
            for i,elem in enumerate(lineelements):
                if i==0:
                    outstring+=elem
                else:
                    outstring+=' '*(spacing[i-1]-len(str(elem)))+str(elem)
            outstring+='\n'
            outfile.write(outstring)
        if chain<nchain-1:
            outfile.write('TER\n')
    outfile.write('END')
    outfile.close()
