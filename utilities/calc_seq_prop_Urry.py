import numpy as np
import matplotlib.pyplot as plt
import os,sys
import scipy.optimize as op

aa_param='''#AA     Mass    Charge  Sigma   Lambda 
ALA     71.08   0.00    5.040   0.602942  
ARG     156.20  1.00    6.560   0.558824
ASN     114.10  0.00    5.680   0.588236
ASP     115.10  -1.00   5.580   0.294119
CYS     103.10  0.00    5.480   0.64706
GLN     128.10  0.00    6.020   0.558824  
GLU     129.10  -1.00   5.920   0.0
GLY     57.05   0.00    4.500   0.57353
HIS     137.10  0.00    6.080   0.764707  
ILE     113.20  0.00    6.180   0.705883
LEU     113.20  0.00    6.180   0.720589
LYS     128.20  1.00    6.360   0.382354
MET     131.20  0.00    6.180   0.676471 
PHE     147.20  0.00    6.360   0.82353 
PRO     97.12   0.00    5.560   0.758824
SER     87.08   0.00    5.180   0.588236 
THR     101.10  0.00    5.620   0.588236 
TRP     186.20  0.00    6.780   0.1 
TYR     163.20  0.00    6.460   0.897059 
VAL     99.07   0.00    5.860   0.664707'''

aa={}
for i in aa_param.split('\n'):
	if i[0]!='#':
		name=i.rsplit()[0]
		other=np.array(i.rsplit()[1:],dtype=float)
		aa[name]=other

seq1to3={'R':'ARG','H':'HIS','K':'LYS','D':'ASP','E':'GLU',
     'S':'SER','T':'THR','N':'ASN','Q':'GLN','C':'CYS',
     'U':'SEC','G':'GLY','P':'PRO','A':'ALA','V':'VAL',
     'I':'ILE','L':'LEU','M':'MET','F':'PHE','Y':'TYR',
     'W':'TRP'}

def seq2para(seq):
	naa=len(list(seq))
	mass=np.zeros(len(seq))
	charge=np.zeros(len(seq))
	sigma=np.zeros(len(seq))
	l=np.zeros(len(seq))
	for idx,i in enumerate(list(seq)):
		i3=seq1to3[i]
		if i3 in aa.keys():
			mass[idx]=aa[i3][0]
			charge[idx]=aa[i3][1]
			sigma[idx]=aa[i3][2]
			l[idx]=aa[i3][3]
		else:
			print(i3)
	return mass,charge,sigma,l
        
def scd(charge):
        res=0.
        for idx,i in enumerate(charge):
                for jdx,j in enumerate(charge):
                        if idx<jdx:
                                res+= i*j*np.sqrt(jdx-idx)
        return res/len(charge)

def shd(l):
        res=0.
        for idx,i in enumerate(l):
                for jdx,j in enumerate(l):
                        if idx<jdx:
                                res+=(i+j)*(jdx-idx)**(-1)
        return res/len(l)
seqfile = open(sys.argv[1]).readlines()
seqlist = []
protnames = []
outfile = open('seq_prop.txt','w')
outfile.write("# Protein list: \n")
for i in seqfile:
    protnames.append(i.split()[0])
    seqlist.append(i.split()[1].strip())
    outfile.write('# %s\n'%i.strip())
outfile.write('# %8s %8s %8s %8s %8s %8s %8s %8s %8s\n'%('SHD(Urry)','SCD','<lambda>(Urry)','net_q','abs_q','f+','f-','length','Mass')) 
for seqno,sequence in enumerate(seqlist):
    #seq = list(open('%s.seq'%sequence).read().strip())
    seq = list(sequence.strip())
    mass,charge,sigma,l=seq2para(seq)
    totmass = np.sum(np.array(mass))
    avel=np.mean(l)
    qp=np.sum(charge[np.where(charge>0)[0]])#/len(charge)
    qm=np.sum(charge[np.where(charge<0)[0]])#/len(charge)
    netq=np.sum(charge)#/len(charge)
    absq=np.sum(np.abs(charge))#/len(charge)

    print(sequence)
    outfile.write('%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f # %s\n'%(shd(l),scd(charge),avel,netq,absq,qp,qm,len(l),totmass,protnames[seqno]))
outfile.close()
