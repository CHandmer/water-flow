#Plot arrays

# User set parameters
zoom=2**4
res = int(45*zoom)
print("res = " + str(res))

# Path to memory location for arrays of a particular resolution
thisdir = "/home/handmer/Documents/Mars/water-flow/"

inputpath = thisdir + "res"+str(res)+"/"

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import image as im
import pandas as pd

place = [2,3]
rawdata = np.load(inputpath+"/test"+str(place[0])+str(place[1])+".npy")[2:-2,2:-2]

mola = rawdata[:,:,0]
depth = rawdata[:,:,2]
#flowNS = np.gradient(mola+depth, axis=0)
flowNS = rawdata[:,:,3]
#flowEW = np.gradient(mola+depth, axis=1)
flowEW = rawdata[:,:,4]
flownorm = rawdata[:,:,7]

totalflow = (flowNS**2+flowEW**2)**0.5


#im.imsave("test.png", flownorm*((0.5+0.5*np.sign(flowNS))*flowNS
#          -(0.5-0.5*np.sign(flowNS))*flowNS))
#im.imsave("test.png", flownorm)
#im.imsave("test.png", totalflow, vmin = 0, vmax = 0.1)
#im.imsave("test.png", depth, vmin =0, vmax = 300)
im.imsave("test.png",np.transpose([mola,depth,totalflow], axes=[1,2,0]))

