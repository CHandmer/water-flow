# This script imports MOLA data from a local directory and creates the arrays at the appropriate resolution.

# This script is modified to increase the ghost zone width from 1 to 2, so that the normalization process works better. It also adds one more layer, to contain the overall delta.

# These arrays are 4+5760/2^n (n = 1, 2, ..., 7) x 9ish, contain
# - MOLA topo data (double)
# - metric (sin latitude)
# - initial depth (initialize at constant GEL, or interpolate/decimate from other file)
# - flow EW
# - flow NS (strictly this is redundant but is a useful piece of data for plotting etc)
# - positive flow, other temp space
# - negative flow
# - norm
# - delta from step to step

# User set parameters
skip = 2**4
res = int(5760/skip)

# Path to memory location for arrays of a particular resolution
thisdir = "/home/handmer/Documents/Mars/water-flow/"
outputpath = thisdir+"res"+str(res)+"/"
sourcepath = thisdir + "RawMola/"
interppath = thisdir + "res180"#"res45-norain-converged"
interpdepth = True

# Data has been pre-sliced, originally derived from
# https://astrogeology.usgs.gov/search/details/Mars/GlobalSurveyor/MOLA/Mars_MGS_MOLA_DEM_mosaic_global_463m/cub
# and sliced into 32 bits
gratextents = [4,8]

# Initialize with global equivalent depth (GED) in meters
GED = 150.

# default order is latitude (N->S), longitude. Watch those square arrays!
#latindex = 0
#lonindex = 0

#from PIL import Image
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

output = np.zeros([res+4,res+4,9])

for latindex in range(gratextents[0]):
    for lonindex in range(gratextents[1]):
        output *= 0.0

        output[2:-2,2:-2,0] = np.loadtxt(sourcepath+"MarsMola"+str(latindex)+"-"+str(lonindex)+".csv", delimiter = ",")[int(skip/2)::skip,int(skip/2)::skip]-2.0**15

        metricvalues = np.sin([np.pi*(i+0.5)/(res*gratextents[0]) for i in range(res*gratextents[0])][latindex*res:(latindex+1)*res])

        for i in range(output.shape[1]-4):
            output[2:-2,i+2,1] = metricvalues

        if interpdepth:
            # Implement some interpolation - constant surface height scheme here. Local tetration.
            interp_input = np.load(interppath+"/test"+str(latindex)+str(lonindex)+".npy")
            output[2:-2:2,2:-2:2,2] = np.maximum(interp_input[2:-2,2:-2,2]+interp_input[2:-2,2:-2,0]-output[2:-2:2,2:-2:2,0],0)
            output[3:-2:2,2:-2:2,2] = np.maximum(interp_input[2:-2,2:-2,2]+interp_input[2:-2,2:-2,0]-output[2:-2:2,2:-2:2,0],0)
            output[2:-2:2,3:-2:2,2] = np.maximum(interp_input[2:-2,2:-2,2]+interp_input[2:-2,2:-2,0]-output[2:-2:2,2:-2:2,0],0)
            output[3:-2:2,3:-2:2,2] = np.maximum(interp_input[2:-2,2:-2,2]+interp_input[2:-2,2:-2,0]-output[2:-2:2,2:-2:2,0],0)
        else:
            output[:,:,2] += GED

        np.save(outputpath+"/test"+str(latindex)+str(lonindex),output)

# reload various files to fill in ghost zone edges
# fix it horizontally first, then vertically in a second pass
for latindex in range(gratextents[0]):
    for lonindex in range(gratextents[1]):
        output = np.load(outputpath+"/test"+str(latindex)+str(lonindex)+".npy")
        output[2:-2,:2] = np.load(outputpath+"/test"+str(latindex)+str((lonindex-1)%gratextents[1])+".npy")[2:-2,-4:-2]
        output[2:-2,-2:] = np.load(outputpath+"/test"+str(latindex)+str((lonindex+1)%gratextents[1])+".npy")[2:-2,2:4]
        np.save(outputpath+"/test"+str(latindex)+str(lonindex),output)
# fix vertically
for latindex in range(gratextents[0]):
    for lonindex in range(gratextents[1]):
        output = np.load(outputpath+"/test"+str(latindex)+str(lonindex)+".npy")
        if latindex == 0:
            # north pole case. replicate previous line, and no water will flow either way.
            output[0] = output[2]
            output[1] = output[2]
        else:
            output[:2] = np.load(outputpath+"/test"+str(latindex-1)+str(lonindex)+".npy")[-4:-2]
        if latindex == gratextents[0]-1:
            #south pole case. replicate polar line, no water will flow either way
            output[-1] = output[-3]
            output[-2] = output[-3]
        else:
            output[-2:] = np.load(outputpath+"/test"+str(latindex+1)+str(lonindex)+".npy")[2:4]
        np.save(outputpath+"/test"+str(latindex)+str(lonindex),output)

