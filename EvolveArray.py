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


import numpy as np

# Create ghost zone arrays. EW, NS, old and new. Contain only depth information.

# Allocate graticule memory

#loop time steps

# refresh ghost zones

# loop over graticules

# load graticule
test = np.load(inputpath+"test00.npy")

# update ghost zones from old ghost zone

# Compute flow for EW in place

# Update graticule depths

# save graticule

# update ghost zones (new)

# Do the same for NS






print(test[:,:,1])
