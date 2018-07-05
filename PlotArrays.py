#Plot arrays

# User set parameters
res = 45

# Path to memory location for arrays of a particular resolution
thisdir = "/home/handmer/Documents/Mars/water-flow/"

inputpath = thisdir + "res"+str(res)+"/"#"-fresh/"

gratextents = [4,8]

import numpy as np
from matplotlib import pyplot as plt

subsample = 1

big_array = np.zeros([int(res*gratextents[0]/subsample),int(res*gratextents[1]/subsample)])

# Allocate graticule memory
# Mola, metric, depth, flowEW, flowNS, positive flow, negative flow, norm
graticule_space = np.zeros([res+2,res+2,6])

# loop over graticules
for latindex in range(gratextents[0]):
    for lonindex in range(gratextents[1]):
        # Load graticule
        graticule_space = np.load(inputpath+"/test"+str(latindex)+str(lonindex)+".npy")
        graticule_space[:,:,5] = (graticule_space[:,:,3]**2+graticule_space[:,:,4]**2)**0.5
        #plt.imshow(graticule_space[1:-1,1:-1,5])
        #plt.show()
        big_array[int(res*latindex/subsample):int(res*(latindex+1)/subsample),
                  int(res*lonindex/subsample):int(res*(lonindex+1)/subsample)] = graticule_space[1:-1,1:-1,2]
        #print([int(res*latindex/subsample),int(res*(latindex+1)/subsample)])

plt.figure(figsize = (20,40))
plt.imshow(big_array)
plt.show()

plt.imshow(np.load(inputpath+"/test"+str(3)+str(0)+".npy")[:,:,0])
plt.show()

