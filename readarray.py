import numpy as np
from PIL import Image
from matplotlib import pyplot as plt


arr = np.load("res7621/array_3_6.npy")

print(arr.shape)

for i in range(3):
    print(np.min(arr[:,:,i]))
    print(np.max(arr[:,:,i]))

#data=arr[2:-2,2:-2,2]
#data = np.interp(data,(data.min(), data.max()), (0,1))

#plt.imshow((0*arr[2:-2,2:-2,0]+arr[2:-2,2:-2,2]))
#plt.show()

#print(data.min())
#print(data.max())

#im = Image.fromarray(data).convert("RGB")
#im.save("testim.png")

#im.show()

deptharray = arr[2:-2,2:-2,2]

depths = np.sort(deptharray.reshape(deptharray.shape[0]*deptharray.shape[1])[::100])

plt.plot(np.log10(depths+10.0**-10))
#plt.show()
plt.savefig("logdepthpowerbig.png")

plt.close("all")
