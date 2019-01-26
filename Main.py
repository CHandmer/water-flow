#This script calls all the relevant python functions. Set these two parameters for the local path and desired global equivalent depth. 

localpath = "./"
# These parameters result in a big northern ocean and lots of lakes, with a number of fluvial canyons entirely flooded. Perhaps half the water is a better starting point?
# A higher degree of precipitation actually leads to more lakes and a smaller ocean.
global_equivalent_depth = 120 #150
precipitation_per_step = 0.0015

# Nothing under here gets changed!

import numpy as np

#Create directory structure.
import subprocess
bashCommand = "mkdir -p res60 res119 res238 res476 res953 res1905 res3811 res7621"
subprocess.check_output(bashCommand, shell=True)

import InitializeArrayFunction
import EvolveArrayFunction

#Initialize res 60
InitializeArrayFunction.InitializeArrays(localpath,60,global_equivalent_depth)

#Evolve res 60 #0.0001,10000
EvolveArrayFunction.EvolveArray(localpath, 60, 0.01, 10000, precipitation_per_step)

# Do the rest of the resolutions
for i in range(1,3,1):#8
    res = int(np.round(7621/2**(7-i)))
    InitializeArrayFunction.InitializeArrays(localpath,res,global_equivalent_depth)
    EvolveArrayFunction.EvolveArray(localpath, res, 0.01, 10000, precipitation_per_step)
