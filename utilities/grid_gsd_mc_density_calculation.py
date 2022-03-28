# Roshan M Regy
# roshanm.regy@tamu.edu
import MDAnalysis as mda 
import numpy as np 
import sys 

# function to fix the pbc
def fixing_pbc(x,l_box):
    l_box_half=l_box/2.0
    x=x-l_box*np.int_(x/l_box_half)
    return x

def get_bins():
    print ("#############Structured binning###########")
    regions = int(input("Enter number of regions: "))
    region_bins = []
    for i in range(regions):
        if i==0:
            start = float(input("STARTing position of region %s: "%(i+1)))
        else:
            start=stop
        stop = float(input("ENDing position of region %s: "%(i+1)))
        nbins = int(input("NUMber of bins in region %s: "%(i+1)))
        region_bins.append(np.linspace(start,stop,nbins+1))
    user_bin_edges = region_bins[0]
    for i,bins in enumerate(region_bins[1:]):
        user_bin_edges = np.concatenate((user_bin_edges,bins[1:]),axis=0)
    return user_bin_edges





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

# Get standard error of mean 
#@jit(nopython=True)
def get_SEM(outmat,rhomat,vol,nblocks,factor):
    bskip=int(rhomat.shape[1]/nblocks)
    for b in np.arange(nblocks):
        outmat[:,b+1] = np.mean(rhomat[:,b*bskip:(b+1)*bskip],axis=1)*factor/vol
    outmat[:,-1] = np.mean(outmat[:,1:nblocks+1],axis=1)
    return outmat

# Process user provided inputs 
def process_input(inputs):
    # Check input arguments 
    if(len(inputs)<9):
        print("All inputs not provided!")
        print ("Inputs: gsd_trajectory indexfile bins start_frame stop_frame skip_frames blocks(SEM) units(nm/ang) choices(if known, otherwise interactive)")  
        exit()

    trajfile = inputs[1]
    univ = mda.Universe(trajfile)
    indexfile = inputs[2]
    bins = int(inputs[3])
    start = int(inputs[4])
    stop = int(inputs[5])
    stride = int(inputs[6])
    noofframes=len(np.arange(start,stop,stride))
    nblocks = int(inputs[7])
    units = inputs[8]
    if units.lower() == 'nm':
        factor = 1.66053904
        outputdist = 1
    elif units.lower() == 'ang':
        factor = 1660.53904
        outputdist = 10
    else:
        print ("Please enter correct units: nm/ang")
        exit()
    # if choices provided with input otherwise ask for choices 
    if len(inputs)>9:
        choices = np.array(inputs[10:]).astype(int)
        index = open(indexfile,'r').readlines()
        indexname = []
        indexlist = []
        for line in index:
            indexname.append(line.split()[0])
            indexlist.append(np.array(line.split()[1:]).astype(int)-1)
        print ("Selected groups for density calculation: ")
        for ch in choices:
            print (indexname[ch])
    elif len(inputs)==9:
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
    return univ,bins,start,stop,stride,noofframes,nblocks,choices,indexlist,indexname,factor,outputdist


### Script execution starts here ###
univ,bins,start,stop,stride,noofframes,nblocks,choices,indexlist,indexname,factor,outputdist = process_input(sys.argv[:])
box_len = univ.dimensions[2]
positions = np.zeros((noofframes,len(univ.atoms.positions[:,-1])))
density = np.zeros((len(choices),bins,noofframes))

# Take all positions from trajectory
for i,ts in enumerate(univ.trajectory[start:stop:stride]):
    positions[i,:] = univ.atoms.positions[:,-1]

# Center the trajectory on maximum density of all atoms
positions = get_centered(positions,bins,noofframes,box_len,univ.atoms[:].masses)

# Calculate density from centered trajectory of each selected group 
user_bin_edges = get_bins()
user_bins = len(user_bin_edges)-1
user_bin_centers = np.zeros(user_bins)
vol = np.zeros(user_bins)
for i,val in enumerate(user_bin_edges[:-1]):
    user_bin_centers[i] = (val+user_bin_edges[i+1])/2.0
    vol[i] = np.absolute(val-user_bin_edges[i+1])*univ.dimensions[0]*univ.dimensions[1]

for ch,choice in enumerate(choices):
    density = np.zeros((user_bins,noofframes))
    outmat = np.zeros((user_bins,nblocks+2))
    for frame in range(noofframes):
        density[:,frame],bin_edges = np.histogram(positions[frame,indexlist[choice]],bins=user_bin_edges,weights=univ.atoms[indexlist[choice]].masses)

    outmat[:,0] = user_bin_centers/outputdist
    outmat = get_SEM(outmat,density,vol,nblocks,factor)
    np.savetxt('%s_z_density.dat'%(indexname[choice]),outmat,fmt='%.3f',header="Bin_center Blocks 1 ... %s Avg_over_blocks"%nblocks)
