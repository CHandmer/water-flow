#This script calls all the relevant python functions. Set these two parameters for the local path and desired global equivalent depth. 

localpath = "./"
# These parameters result in a big northern ocean and lots of lakes, with a number of fluvial canyons entirely flooded. Perhaps half the water is a better starting point?
# A higher degree of precipitation actually leads to more lakes and a smaller ocean.
global_equivalent_depth = 120 #150
precipitation_per_step = 0.0015

# Nothing under here gets changed!


#Create directory structure.
import subprocess
bashCommand = "mkdir -p res45 res90 res180 res360 res720 res1440 res2880 res5760"
subprocess.check_output(bashCommand, shell=True)

import InitializeArrayFunction
import EvolveArrayFunction

#Initialize res 45
InitializeArrayFunction.InitializeArrays(localpath,45,global_equivalent_depth)

#Evolve res 45
EvolveArrayFunction.EvolveArray(localpath, 45, 0.0001, 10000, precipitation_per_step)

# Do the rest of the resolutions
for i in range(1,8,1):
    res = 45*2**i
    InitializeArrayFunction.InitializeArrays(localpath,res,global_equivalent_depth)
    EvolveArrayFunction.EvolveArray(localpath, res, 0.001, 1000, precipitation_per_step)
