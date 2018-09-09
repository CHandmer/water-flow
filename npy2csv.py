import numpy as np

inputpath = "/home/handmer/Documents/Mars/water-flow/"

rawdata = np.load(inputpath+"/gale_"+str(30.0)+"_"+str(0.0)+"_"+str(0)+".npy")

np.savetxt(inputpath+"/depth4.csv", rawdata[2:-2:4,2:-2:4,2], delimiter=",")
