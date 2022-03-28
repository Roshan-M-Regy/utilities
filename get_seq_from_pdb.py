# Roshan M Regy
# roshanm.regy@tamu.edu
import sys

amino_3to1 = {
        'ALA':'A',
        'CYS':'C',
        'ASP':'D',
        'GLU':'E',
        'PHE':'F',
        'GLY':'G', 
        'HIS':'H', 
        'ILE':'I', 
        'LYS':'K', 
        'LEU':'L', 
        'MET':'M', 
        'ASN':'N', 
        'PRO':'P',  
        'GLN':'Q', 
        'ARG':'R', 
        'SER':'S', 
        'THR':'T', 
        'VAL':'V',
        'TRP':'W',
        'TYR':'Y', 
        }
lines = open(sys.argv[1]).readlines()
seq = ''
seqfile = open('pdbseq.txt','w')
for line in lines:
    if len(line.split())==0:
        continue
    elif line.split()[0] not in ['ATOM','TER','END']:
        continue
    elif line.split()[0] == 'TER':
        seq += '\n'
        seqfile.write(seq)
        seq = ''
    elif line.split()[0] == 'END':
        if seq!='':
            seq += '\n'
            seqfile.write(seq)
            seq = ''
        break
    else:
            if line.split()[2]=='CA':
                seq+=amino_3to1[line.split()[3]]
seqfile.close()
