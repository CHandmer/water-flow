# This script contains functions of the primitive array that evolve it forward one time step. 

# Steps
# 1) Read array from memory
# 2) Update ghost zones for outer edges
# 3) Evolve array (EW)
# 4) Update ghost zones for inner edges to new ghost zone array
# 5) Write update to memory

# 5.5) set new ghost zones to old ghost zones
# 6) Do it all again for NS.


# This script also initializes an array containing ghost zone information, which is how adjacent sectors pass depth information back and forth.

# User set parameters
res = 45

# Path to memory location for arrays of a particular resolution
thisdir = "/home/handmer/Documents/Mars/water-flow/"

inputpath = thisdir + "res"+str(res)+"/"

gratextents = [4,8]

timestep = 0.1

precip = 0*0.0015*gratextents[0]*res/720

reset_depths = False
GED = 150

import numpy as np
from matplotlib import pyplot as plt

# Create ghost zone arrays. EW, NS, old and new. Contain only depth information.
# Each graticule has a one element ghost zone. So the ghost zone for NS work must be 2*gratextents[0] x res*gratextents[1] and EW must be gratextents[0]*res x 2*gratextents[1]
# It doesn't need to include the diagonals because they're irrelevant for the algo.
# Note to self: This would have been much easier to implement as methods on a class. Live and learn. 
ghostNSold = np.zeros([2*gratextents[0],res*gratextents[1]])
ghostNSnew = np.zeros([2*gratextents[0],res*gratextents[1]])
ghostEWold = np.zeros([res*gratextents[0],2*gratextents[1]])
ghostEWnew = np.zeros([res*gratextents[0],2*gratextents[1]])

# Allocate graticule memory
# Mola, metric, depth, flowEW, flowNS, positive flow, negative flow, norm
graticule_space = np.zeros([res+2,res+2,6])

# Fill in ghost zone depth numbers (initialization)
# Imagine there are two domains, one next to the other, and not periodic boundaries (zero flow)
# 1* 1 1 1 1 g1
#          g2 2 2 2 2 2*
# When domain 1 runs, it duplicates its first row to make 1*, then reads the first row of 2 to make g1. 
# When it finishes, it writes its rightmost value to g2. So there's a "crossing of arms" between reading and writing. 
# The only memory structure that will make any sense in the long term is one that is only ghost zones in the correct order.

# For ghostNS, it only needs 6 rows, since the poles are redundant, but 8 rows simplifies addressing.
# The 0th and 1st row are the interior 1 and -2th rows of domain latindex = 0, and so on .

# Image four domains, next to each other, periodic boundaries
# g1 1 1 1 1 g1     g3 3 3 3 3 g3     g1 1 1 1 1 g1
#          g2 2 2 2 2 g2     g4 4 4 4 4 g4
# Domain 1 reads 4 and 2 to get g1s, which requires some cyclic references. If the ghost zones are labeled in order, they are g4, g2, g1, g3, g2, g4, g3, g1, ...

# In each cycle, it will read remote ghost zones, then write local ones. It's a beautiful thing. 

# minalt allows the precipitation algorithm to normalize for heights. 
# 
minalt = 0
totalalts = 0
for latindex in range(gratextents[0]):
    for lonindex in range(gratextents[1]):
        graticule_space = np.load(inputpath+"/test"+str(latindex)+str(lonindex)+".npy")
        if reset_depths:
            graticule_space[:,:,2] *= 0
            graticule_space[:,:,2] += GED
        ghostNSnew[2*latindex,res*lonindex:res*(lonindex+1)] = graticule_space[1,1:-1,2]
        ghostNSnew[1+2*latindex,res*lonindex:res*(lonindex+1)] = graticule_space[-2,1:-1,2]
        ghostEWnew[res*latindex:res*(latindex+1),2*lonindex] = graticule_space[1:-1,1,2]
        ghostEWnew[res*latindex:res*(latindex+1),1+2*lonindex] = graticule_space[1:-1,-2,2]

        # set minalt
        minalt = np.min([np.min(graticule_space[:,:,0]),minalt])
        totalalts += np.sum(graticule_space[1:-1,1:-1,0])

# rescale to set the minimum altitude to "0"
totalalts += -res*res*gratextents[0]*gratextents[1]*minalt

# This is an array of numbers representing the total flow within the system. As a measure of convergence.
total_flow = []

#loop time steps
for i in range(1):

    # refresh ghost zones (begin the new time step with fresh space for old ghost zones)
    ghostEWold = ghostEWnew

    # refresh total evaporation to zero
    total_evap = 0

    # add a new term to total flow. At this point, doing EW and NS flows separately. 
    total_flow.append(0)

    # loop over graticules
    for latindex in range(gratextents[0]):
        for lonindex in range(gratextents[1]):
            # Load graticule
            graticule_space = np.load(inputpath+"/test"+str(latindex)+str(lonindex)+".npy")
            
            # update ghost zones from old ghost zone 
            # (which was the new ghost zone from the previous time step)
            graticule_space[1:-1,0,2] = ghostEWold[res*latindex:res*(latindex+1),(2*lonindex-1)%gratextents[1]]
            graticule_space[1:-1,-1,2] = ghostEWold[res*latindex:res*(latindex+1),(2*lonindex+2)%gratextents[1]]
            
            # Compute flow for EW in place
            graticule_space[:,:-1,3] = timestep*np.diff(graticule_space[:,:,0] + graticule_space[:,:,2], axis = 1)
            
            # Update graticule depths
            # This must be normalized by the existing depth, or depth will go negative. 
            # Important to preserve conservatism.
            # Compute positive flow (right to left)
            graticule_space[:,:-1,5] = (0.5+0.5*np.sign(graticule_space[:,:-1,3]))*graticule_space[:,:-1,3] #positive, flows to the left, loss from right, add to left
            # Compute negative flow (left to right)
            graticule_space[:,:-1,6] = (0.5-0.5*np.sign(graticule_space[:,:-1,3]))*graticule_space[:,:-1,3] #negative, flows to the right, loss from left, add to right
            # Compute norm
            graticule_space[:,1:-1,7] = np.minimum(graticule_space[:,1:-1,2]/(10**-8 + graticule_space[:,:-2,5] - graticule_space[:,1:-1,6]),1.0)
            
            # update depth, including a metric term that converts flows to volumes, scaled by latitude.
            graticule_space[:,1:-1,2] += graticule_space[:,1:-1,7]*(graticule_space[:,1:-1,5] - graticule_space[:,:-2,5] + graticule_space[:,1:-1,6] - graticule_space[:,:-2,6])*graticule_space[:,1:-1,1]
            total_flow[-1] += np.sum(np.abs(graticule_space[:,1:-1,7]*(graticule_space[:,1:-1,5] - graticule_space[:,:-2,5] + graticule_space[:,1:-1,6] - graticule_space[:,:-2,6])*graticule_space[:,1:-1,1]))
            
            # Note to self. Imagine an array like 
            # 0 0 1 0 0 0
            # Diff looks like
            # 0 1 -1 0 0 
            # Compute norm fraction for central 4
            # Diff+[:-1] - Diff-[1:]
            # Compute fate for central 4
            # (Diff+[1:] - Diff+[:-1] + Diff-[1:] - Diff-[:-1]) * norm
            
            # Determine correct evaporation level. Reusing layer 5 of graticule array.
            graticule_space[1:-1,1:-1,5] = np.minimum(graticule_space[1:-1,1:-1,2],precip)
            # Evaporate
            graticule_space[1:-1,1:-1,2] -= graticule_space[1:-1,1:-1,5]
            # Compute total volume by including metric
            total_evap += np.sum(graticule_space[1:-1,1:-1,5]*graticule_space[1:-1,1:-1,1])
            
            # update ghost zones (new)
            ghostEWnew[res*latindex:res*(latindex+1),2*lonindex] = graticule_space[1:-1,1,2]
            ghostEWnew[res*latindex:res*(latindex+1),1+2*lonindex] = graticule_space[1:-1,-2,2]
            
            # save graticule
            np.save(inputpath+"/test"+str(latindex)+str(lonindex),graticule_space)



    # Do the same for NS
    # Cycle ghostzone. I hope this is assignation, not binding.
    ghostNSold = ghostNSnew
    
    total_flow.append(0)
    
    # loop over graticules
    for latindex in range(gratextents[0]):
        for lonindex in range(gratextents[1]):
            # Load graticule
            graticule_space = np.load(inputpath+"/test"+str(latindex)+str(lonindex)+".npy")
            
            # update ghost zones from old ghost zone 
            # (which was the new ghost zone from the previous time step)
            if latindex == 0:
                graticule_space[0,1:-1,2] = graticule_space[1,1:-1,2]
            else:
                graticule_space[0,1:-1,2] = ghostNSold[2*latindex-1,res*lonindex:res*(lonindex+1)]
            if latindex == gratextents[0]-1:
                graticule_space[-1,1:-1,2] = graticule_space[-2,1:-1,2]
            else:
                graticule_space[-1,1:-1,2] = ghostNSold[2*latindex+2,res*lonindex:res*(lonindex+1)]
        
            # Rain including metric
            graticule_space[1:-1,1:-1,2] += total_evap*graticule_space[1:-1,1:-1,0]/(totalalts*graticule_space[1:-1,1:-1,1])

            # Compute flow for NS in place
            graticule_space[:-1,:,4] = timestep*np.diff(graticule_space[:,:,0] + graticule_space[:,:,2], axis = 0)

            # Update graticule depths
            # This must be normalized by the existing depth, or depth will go negative. 
            # Important to preserve conservatism.
            # Compute positive flow (right to left)
            graticule_space[:-1,:,5] = (0.5+0.5*np.sign(graticule_space[:-1,:,4]))*graticule_space[:-1,:,4] #positive, flows to the left, loss from right, add to left
            # Compute negative flow (left to right)
            graticule_space[:-1,:,6] = (0.5-0.5*np.sign(graticule_space[:-1,:,4]))*graticule_space[:-1,:,4] #negative, flows to the right, loss from left, add to right
            # Compute norm
            graticule_space[1:-1,:,7] = np.minimum(graticule_space[1:-1,:,2]/(10**-8 + graticule_space[:-2,:,5] - graticule_space[1:-1,:,6]),1.0)

            # update depth, including a metric term that converts flows to volumes, scaled by latitude.
            graticule_space[1:-1,:,2] += graticule_space[1:-1,:,7]*(graticule_space[1:-1,:,5] - graticule_space[:-2,:,5] + graticule_space[1:-1,:,6] - graticule_space[:-2,:,6])*graticule_space[1:-1,:,1]
            total_flow[-1] += np.sum(np.abs(graticule_space[1:-1,:,7]*(graticule_space[1:-1,:,5] - graticule_space[:-2,:,5] + graticule_space[1:-1,:,6] - graticule_space[:-2,:,6])*graticule_space[1:-1,:,1]))

            # update ghost zones (new)
            ghostNSnew[2*latindex,res*lonindex:res*(lonindex+1)] = graticule_space[1,1:-1,2]
            ghostNSnew[1+2*latindex,res*lonindex:res*(lonindex+1)] = graticule_space[-2,1:-1,2]

            # save graticule
            np.save(inputpath+"/test"+str(latindex)+str(lonindex),graticule_space)
            
#plt.plot((np.array(total_flow)[::2]**2+np.array(total_flow)[1::2]**2)**0.5)
#plt.show()

