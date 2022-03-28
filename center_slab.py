import MDAnalysis as mda 
import numpy as np 
import sys 
import gsd.hoomd
# function to fix the pbc
def fixing_pbc(x,l_box):
    l_box_half=l_box/2.0
    x=x-l_box*np.int_(x/l_box_half)
    return x
# Center frames 
def get_centered(positions,n_bins,noofframes,box_len,masses):
    bins = np.linspace(0,box_len,n_bins+1)
    bin_w = bins[1]-bins[0]
    density, bin_edges = np.histogram(positions[0,:]+box_len/2.0,bins,weights=masses)
    xavg = density.mean()
    for frame in (range(noofframes)):
            #box size and area in each dimension
            zvalues = positions[frame,:]+box_len/2. # shif positions so they start from zero and end at box_len
            density, bin_edges = np.histogram(zvalues,bins,weights=masses)
            delxA = density - xavg
            ## Get autocorrelation
            C_AA = np.zeros(n_bins)
            C_AA[0] = (delxA**2).sum()
            for d in range(1, n_bins):
                C_AA[d] = (delxA*np.concatenate([delxA[d:], delxA[:d]], axis=0)).sum()
            #ignors confs that lead to C_AA=0
            if((C_AA.max() - C_AA.min())==0.0):
                print(density)
                continue
                print("test", C_AA.max())
            C_AA /= ((C_AA.max() - C_AA.min())/2.)
            C_AA -= (C_AA.max() - 1) # Normalize to range from -1 to 1
            C_AA = np.concatenate([C_AA[int(n_bins/2):],C_AA[:int(n_bins/2)]], axis=0)
            ## Get pulse function
            p_A = C_AA > 0

            ## Get overlap function
            X_c = np.zeros(n_bins)
            X_c[0] = (density*p_A).sum()
            for n in range(1, n_bins):
                X_c[n] = (density*np.concatenate([p_A[n:], p_A[:n]], axis=0)).sum()
            X_c = np.concatenate([X_c[int(n_bins/2):], X_c[:int(n_bins/2)]], axis=0)

            ## Get number of points where maximum, and select one in center
            count = (X_c == X_c.max()).sum()
            n_max = (X_c == X_c.max()).nonzero()[0][int(count/2)]
            positions[frame] = positions[frame]+(n_max)*bin_w-box_len/2.0
            positions[frame] = fixing_pbc(positions[frame],box_len)
    
    return positions

# Process user provided inputs 
def process_input(inputs):
    # Check input arguments 
    print ("len of inputs")
    print (len(inputs))
    if(len(inputs)<6):
        print("All inputs not provided!")
        print ("Inputs: dcd/xtc_trajectory topfile(.data) indexfile bins units")  
        exit()

    trajfile = inputs[1]
    topfile = inputs[2]
    univ = mda.Universe(topfile,trajfile)
    indexfile = inputs[3]
    bins = int(inputs[4])
    noofframes=len(univ.trajectory)
    units = inputs[5]
    if units.lower() == 'nm':
        factor = 1.66053904*1000
        outputdist = 1#10
    elif units.lower() == 'ang':
        factor = 1660.53904/1000
        outputdist = 1#1
    else:
        print ("Please enter correct units: nm/ang")
        exit()
    # if choices provided with input otherwise ask for choices 
    if len(inputs)>6:
        choices = np.array(inputs[6:]).astype(int)
        index = open(indexfile,'r').readlines()
        indexname = []
        indexlist = []
        for line in index:
            indexname.append(line.split()[0])
            indexlist.append(np.array(line.split()[1:]).astype(int)-1)
        print ("Selected groups for density calculation: ")
        for ch in choices:
            print (indexname[ch])
    elif len(inputs)==6:
        index = open(indexfile,'r').readlines()
        indexname = ['All']
        indexlist = [range(len(univ.atoms.positions[:,1]))]
        print ("Index")
        print ("0: All")
        for h,line in enumerate(index):
            indexname.append(line.split()[0])
            print ("%s: %s"%(h+1,indexname[-1]))
            indexlist.append(np.array(line.split()[1:]).astype(int)-1)
        choices = np.array(input("Enter choices for density calculation: ").split()).astype(int)
    return univ,bins,noofframes,choices,indexlist,indexname,factor,outputdist
### Script execution starts here ###
univ,bins,noofframes,choices,indexlist,indexname,factor,outputdist = process_input(sys.argv[:])

positions = np.zeros((noofframes, univ.atoms.positions.shape[0]))
allcoords = np.zeros((noofframes, univ.atoms.positions.shape[0],univ.atoms.positions.shape[1]))

for dim in range(3):
    # Take all positions from trajectory
    for i,ts in enumerate(univ.trajectory):
        positions[i,:] = univ.atoms.positions[:,dim]
    box_len = univ.dimensions[dim]
    # Center the trajectory on maximum density of all atoms
    if dim==2:
        positions = get_centered(positions,bins,noofframes,box_len,univ.atoms[:].masses)
    allcoords[:,:,dim] = positions
'''
newgsd = gsd.hoomd.open('testtraj.gsd','wb')
oldgsd = gsd.hoomd.open('prod.gsd','rb')[0]
for i in range(noofframes):
    oldgsd.particles.position = allcoords[i,:,:]*outputdist
    newgsd.append(oldgsd)
newgsd.close()
'''
w = mda.Writer(sys.argv[1][:-4]+'_centered.xtc',univ.atoms)
w.dt = 1000
for i,ts in enumerate(univ.trajectory):
    frame = ts
    frame.positions = allcoords[i,:,:]*outputdist
    frame._unitcell *= outputdist
    w.write(univ.atoms)
w.close()
