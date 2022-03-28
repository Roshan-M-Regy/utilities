# Roshan M Regy
# roshanm.regy@tamu.edu
import sys

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

lines = open(sys.argv[1]).readlines()
seq = open(sys.argv[2]).readlines()

newpdbfile  = open('newpdb.pdb','w')
chaincounter = 0
aacounter = 0
for i,line in enumerate(lines):
    if len(line.split())<1:
        continue
    elif line.split()[0]=='TER':
        newpdbfile.write(line)
        chaincounter+=1
        aacounter = 0
    elif line.split()[0]=='END':
        newpdbfile.write(line)
        break

    else:
        if line.split()[2]=='CA':
            print (line)
            for ch,charc in enumerate(line):
                if line[ch:ch+3] in amino_1to3.values():
                    print (line[ch:ch+3])
                    break
            line=list(line)
            line[ch:ch+3] = str(amino_1to3[seq[chaincounter][aacounter]])
            newpdbfile.write(''.join(line))
            aacounter+=1

newpdbfile.close()
        
