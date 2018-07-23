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
import colorsys as cs

place = [2,3]
rawdata = np.load(inputpath+"/array"+str(place[0])+str(place[1])+".npy")[2:-2,2:-2]

mola = rawdata[:,:,0]
depth = rawdata[:,:,2]
flowNS2 = np.gradient(mola+depth, axis=0)
flowNS = rawdata[:,:,3]
flowEW2 = np.gradient(mola+depth, axis=1)
flowEW = rawdata[:,:,4]
flownorm = rawdata[:,:,7]

totalflow = (flowNS**2+flowEW**2)**0.5
totalflownorm = totalflow/np.max(totalflow)
flowdir = np.arctan2(flowNS2,flowEW2)/(2*np.pi)

relief = np.gradient(mola, axis=0) + np.gradient(mola, axis=1)
reliefnorm = relief/np.max(relief)

relief2 = np.gradient(mola+depth, axis=0) + np.gradient(mola+depth, axis=1)
relief2norm = relief2/np.max(relief2)

steepness = (np.gradient(mola, axis=0)**2 + np.gradient(mola, axis=1)**2)**0.5
steepnessnorm = steepness/np.max(steepness)

#im.imsave("test.png", flownorm*((0.5+0.5*np.sign(flowNS))*flowNS
#          -(0.5-0.5*np.sign(flowNS))*flowNS))
#im.imsave("test.png", flownorm)
#im.imsave("test.png", totalflow, vmin = 0, vmax = 0.1)
#im.imsave("test.png", depth, vmin =0, vmax = 300)
#im.imsave("test.png",np.transpose([mola,depth,totalflow], axes=[1,2,0]))
#im.imsave("test.png",relief)

#im.imsave("test.png",cs.hsv_to_rgb(flowdir.all(), totalflownorm.all(), reliefnorm.all()))

#im.imsave("test.png",np.array([[cs.hsv_to_rgb(flowdir[i,j], 0.1+10*totalflownorm[i,j], 0.5+0.4*reliefnorm[i,j]) for i in range(flowdir.shape[0])] for j in range(flowdir.shape[1])]))

#im.imsave("test.png",np.array([[cs.hsv_to_rgb(steepnessnorm[i,j]**0.25,
#                                              0.5+0.5*np.sign(relief2norm[i,j])*(np.sign(relief2norm[i,j])*relief2norm[i,j])**0.5,
#                                              totalflownorm[i,j]**0.25
#                                          ) for i in range(flowdir.shape[0])] for j in range(flowdir.shape[1])]))


# In order to do this semi-sensibly, I need to know the distribution of these various properties. I also need a better colormap.

#plt.plot(np.sort(steepnessnorm[::10,::10].reshape(72*72))**0.25)
#plt.show()

#rnsign = np.sign(reliefnorm[::10,::10])
#plt.plot(np.sort((rnsign*(rnsign*reliefnorm[::10,::10])**1).reshape(72*72)))
#plt.show()

#plt.plot(np.tanh(10*np.sort(reliefnorm[::10,::10].reshape(72*72))))
#plt.show()

#plt.plot(np.tanh(100*np.sort(totalflownorm[::10,::10].reshape(72*72))))
#plt.show()

#plt.plot(np.tanh(0.002*np.sort(depth[::10,::10].reshape(72*72))))
#plt.show()


#im.imsave("test.png",np.array([[cs.hsv_to_rgb(1-0.3*np.tanh(0.002*depth[i,j]),
#                                              0.5-0.4*np.tanh(10*reliefnorm[i,j]),
#                                              0.4+0.6*np.tanh(100*totalflownorm[i,j])
#                                          ) for i in range(flowdir.shape[0])] for j in range(flowdir.shape[1])]))

#im.imsave("test.png",np.array([[cs.hsv_to_rgb(1-0.6*np.tanh(0.002*depth[i,j]),
#                                              0.4+0.6*np.tanh(100*totalflownorm[i,j]),
#                                              0.5-0.4*np.tanh(10*reliefnorm[i,j])
#                                          ) for i in range(flowdir.shape[0])] for j in range(flowdir.shape[1])]))

#Time to home cook a HSB2RGB function.
# from https://peterkovesi.com/projects/colourmaps/
from PIL import Image
colorfunctiondata = Image.open("diverging-linear_bjr_30-55_c53_n256.png")
cfdata = np.asarray(colorfunctiondata)[-1]
#print(cfdata.shape)

def MyHSB2RGBfunc(hue, sat, bri):
    #takes and outputs values between 0 and 1
    hue = cfdata[int(round(511*hue))]/255.
    hueav = np.sum(hue)/3
    totalsat = np.array([hueav for i in range(3)])
    return bri*(sat*hue + (1-sat)*totalsat)

#print([[0.1*depthindex,round(511*np.tanh(0.002*0.1*depthindex))] for depthindex in range(1001)])

#im.imsave("test.png",np.array([[MyHSB2RGBfunc(1-np.tanh(0.1*0.1*depthindex),0.01*satindex,1) for depthindex in range(501)] for satindex in range(101)]))



#im.imsave("test.png",np.array([[MyHSB2RGBfunc(1-np.tanh(0.01*depth[i,j]),
#                                              0.2+0.8*np.tanh(100*totalflownorm[i,j]),
#                                              0.6+0.4*np.tanh(10*reliefnorm[i,j])
#                                          ) for i in range(flowdir.shape[0])] for j in range(flowdir.shape[1])]))

#im.imsave("test.png",np.array([[MyHSB2RGBfunc(1-0*np.tanh(0.1*depth[i,j]),
#                                              0.5-0.5*np.tanh(10*reliefnorm[i,j]),
#                                              0.5+0.5*np.tanh(100*totalflownorm[i,j])
#                                          ) for i in range(flowdir.shape[0])] for j in range(flowdir.shape[1])]))


#im.imsave("test.png",np.array([[MyHSB2RGBfunc(1-np.tanh(0.01*depth[i,j]),
#                                              0.5+0.25*(np.tanh(100*totalflownorm[i,j])-np.tanh(10*reliefnorm[i,j])),
#                                              0.5+0.25*(np.tanh(10*reliefnorm[i,j])+np.tanh(100*totalflownorm[i,j]))
#                                          ) for i in range(flowdir.shape[0])] for j in range(flowdir.shape[1])]))

#This is a good one
im.imsave("test.png",
          np.array([[MyHSB2RGBfunc(0.9-0.7*np.tanh(0.01*depth[i,j]),
                                   1-np.tanh(100*totalflownorm[i,j]),
                                   0.5+0.4*np.tanh(np.tanh(10*reliefnorm[i,j])+np.tanh(100*totalflownorm[i,j]))
                               ) 
                     for i in range(flowdir.shape[0])] 
                    for j in range(flowdir.shape[1])]))


# Let's try for flow direction too.
colorfunctiondata2 = Image.open("cyclic-isoluminant.png")
cfdata2 = np.asarray(colorfunctiondata2)[-3,10:-5]
#print(cfdata2.shape)

def MyHSB2RGBfunc2(hue, sat, bri):
    #takes and outputs values between 0 and 1
    hue = cfdata2[int(round((cfdata2.shape[0]-1)*hue))]/255.
    hueav = np.sum(hue)/3
    totalsat = np.array([hueav for i in range(3)])
    return bri*(sat*hue + (1-sat)*totalsat)


#im.imsave("test.png",
#          np.array([[MyHSB2RGBfunc(0.9-0.7*np.tanh(0.01*depth[i,j]),
#                                   1-np.tanh(100*totalflownorm[i,j]),
#                                   0.5+0.4*np.tanh(np.tanh(10*reliefnorm[i,j])+np.tanh(100*totalflownorm[i,j])))*(1-np.tanh(0.1*depth[i,j])) + 
#                     MyHSB2RGBfunc2(0.5+flowdir[i,j],
#                                    1-np.tanh(50*totalflownorm[i,j]),
#                                    0.5+0.4*np.tanh(np.tanh(10*reliefnorm[i,j])+np.tanh(100*totalflownorm[i,j])))*np.tanh(0.1*depth[i,j])
#                     for i in range(flowdir.shape[0])] 
#                    for j in range(flowdir.shape[1])]))
