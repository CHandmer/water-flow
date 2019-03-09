#This script calls all the relevant python functions. Set these two parameters for the local path and desired global equivalent depth.

#Before running this run the wget, splittif and readtif scripts to put the giant DEM file in the correct form.

localpath = "./"
# These parameters result in a big northern ocean and lots of lakes, with a number of fluvial canyons entirely flooded. Perhaps half the water is a better starting point?
# A higher degree of precipitation actually leads to more lakes and a smaller ocean.
global_equivalent_depth = 50 #120 #150 # Perhaps 50m is closer to the truth?
# I want to adjust this so peak precip is 1/50th of 0.0015, and scales with resolution.
precipitation_per_step = 0.0015 # This might be a bit high, not sure yet.
#precipitation_per_step = 0.004 # This might be a bit high, not sure yet.
# Given that more steps are required to flow a given distance at higher resolution, a conservative approach would be to reduce precipitation_per_step as resolution increases. I have decided not to because, given rather coarse resolution, increasing rainfall as time goes on tends to fully charge and render visible the tiny tributaries, making the map more interesting.

# Nothing under here gets changed! Hahaha.

import numpy as np

#Create directory structure.
import subprocess
bashCommand = "mkdir -p res60 res119 res238 res476 res953 res1905 res3810 res7621"
subprocess.check_output(bashCommand, shell=True)

import InitializeArrayFunction
import EvolveArrayFunction

#Initialize res 60
InitializeArrayFunction.InitializeArrays(localpath,60,global_equivalent_depth)

#Evolve res 60 #0.0001,10000. This may be overly conservative. 
EvolveArrayFunction.EvolveArray(localpath, 60, 0.0001, 10000, precipitation_per_step)

# Do the rest of the resolutions
maxsteps=10000
for i in range(1,8,1):
    res = int(np.round(7621/2**(7-i)))
    InitializeArrayFunction.InitializeArrays(localpath,res,global_equivalent_depth)
    if i==7:
        maxsteps=10 # This prevents all the water disappearing into gullies in the last step
    EvolveArrayFunction.EvolveArray(localpath, res, 0.001, maxsteps, precipitation_per_step)#/2**i)
    # Divide precip by 2 every subdivision to capture the effective timestep halving. 2^7 = 128, so final drainage should be much better. Keep an eye on it. 
