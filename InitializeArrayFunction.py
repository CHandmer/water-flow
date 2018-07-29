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

#from PIL import Image
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import pandas as pd
from matplotlib import pyplot as plt

def InitializeArrays(path,resolution,initialdepth):
    res = resolution
    print("res = " + str(res))

    skip = int(5760/res)

    # Path to memory location for arrays of a particular resolution
    thisdir = path
    outputpath = thisdir+"res"+str(res)+"/"
    sourcepath = thisdir + "RawMola/"
    interppath = thisdir + "res"+str(int(res/2))#"res2880"#"res1440"#"res720"#"res360"#"res180"#"res90"#"res45-norain-converged"
    if res==45:
        interpdepth = False
    else:
        interpdepth = True

    # Data has been pre-sliced, originally derived from
    # https://astrogeology.usgs.gov/search/details/Mars/GlobalSurveyor/MOLA/Mars_MGS_MOLA_DEM_mosaic_global_463m/cub
    # and sliced into 32 bits
    gratextents = [4,8]

    # Initialize with global equivalent depth (GED) in meters
    GED = initialdepth

    # default order is latitude (N->S), longitude. Watch those square arrays!

    output = np.zeros([res+4,res+4,9])

    for latindex in range(gratextents[0]):
        for lonindex in range(gratextents[1]):
            print("Initializing array ("+str(latindex)+","+str(lonindex)+")")

            output *= 0.0

            output[2:-2,2:-2,0] = np.loadtxt(sourcepath+"MarsMola"+str(latindex)+"-"+str(lonindex)+".csv", delimiter = ",")[int(skip/2)::skip,int(skip/2)::skip]-2.0**15

            metricvalues = np.sin([np.pi*(i+0.5)/(res*gratextents[0]) for i in range(res*gratextents[0])][latindex*res:(latindex+1)*res])

            for i in range(output.shape[1]-4):
                output[2:-2,i+2,1] = metricvalues

            if interpdepth:
                # This interpolation scheme, while creating slightly bumpy ocean, avoids the creation of huge quantities of unwanted water on dry slopes
                interp_input = np.load(interppath+"/array"+str(latindex)+str(lonindex)+".npy")#This generates unwanted water on steep, dry, slopes. 
                output[2:-2:2,2:-2:2,2] = np.maximum(interp_input[2:-2,2:-2,2],0)#+interp_input[2:-2,2:-2,0]-output[2:-2:2,2:-2:2,0],0)
                output[3:-2:2,2:-2:2,2] = np.maximum(interp_input[2:-2,2:-2,2],0)#+interp_input[2:-2,2:-2,0]-output[3:-2:2,2:-2:2,0],0)
                output[2:-2:2,3:-2:2,2] = np.maximum(interp_input[2:-2,2:-2,2],0)#+interp_input[2:-2,2:-2,0]-output[2:-2:2,3:-2:2,0],0)
                output[3:-2:2,3:-2:2,2] = np.maximum(interp_input[2:-2,2:-2,2],0)#+interp_input[2:-2,2:-2,0]-output[3:-2:2,3:-2:2,0],0)
            else:
                output[:,:,2] += GED
                
            np.save(outputpath+"/array"+str(latindex)+str(lonindex),output)

    # reload various files to fill in ghost zone edges
    # fix it horizontally first, then vertically in a second pass
    print("Completing ghost zones")
    for latindex in range(gratextents[0]):
        for lonindex in range(gratextents[1]):
            output = np.load(outputpath+"/array"+str(latindex)+str(lonindex)+".npy")
            output[2:-2,:2] = np.load(outputpath+"/array"+str(latindex)+str((lonindex-1)%gratextents[1])+".npy")[2:-2,-4:-2]
            output[2:-2,-2:] = np.load(outputpath+"/array"+str(latindex)+str((lonindex+1)%gratextents[1])+".npy")[2:-2,2:4]
            np.save(outputpath+"/array"+str(latindex)+str(lonindex),output)
    # fix vertically
    for latindex in range(gratextents[0]):
        for lonindex in range(gratextents[1]):
            output = np.load(outputpath+"/array"+str(latindex)+str(lonindex)+".npy")
            if latindex == 0:
                # north pole case. replicate previous line, and no water will flow either way.
                output[0] = output[2]
                output[1] = output[2]
            else:
                output[:2] = np.load(outputpath+"/array"+str(latindex-1)+str(lonindex)+".npy")[-4:-2]
            if latindex == gratextents[0]-1:
                #south pole case. replicate polar line, no water will flow either way
                output[-1] = output[-3]
                output[-2] = output[-3]
            else:
                output[-2:] = np.load(outputpath+"/array"+str(latindex+1)+str(lonindex)+".npy")[2:4]
            np.save(outputpath+"/array"+str(latindex)+str(lonindex),output)
    print("Initialization complete")

