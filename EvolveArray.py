# This script contains functions of the primitive array that evolve it forward one time step. 

# Steps
# 1) Read array from memory
# 2) Update ghost zones for outer edges
# 3) Evolve array (EW)
# 4) Update ghost zones for inner edges to new ghost zone array
# 5) Write update to memory

# 5.5) set new ghost zones to old ghost zones
# 6) Do it all again for NS.
