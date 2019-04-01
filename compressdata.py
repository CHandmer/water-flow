import numpy as np
import os

os.system("mkdir res7621_compressed")

def compressArray(lati,loni):
    data = np.load("res7621/array_"+str(lati)+"_"+str(loni)+".npy")
    print("saving "+str(lati)+"_"+str(loni))
    np.save("res7621_compressed/array_"+str(lati)+"_"+str(loni)+".npy", data[:,:,0:3:2])

for lati in range(7):
    for loni in range(14):
        compressArray(lati,loni)

    
 
