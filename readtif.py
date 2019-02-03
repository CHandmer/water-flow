from PIL import Image
import numpy as np

for i in range(14):
    for j in range(7):
        np.save("output_dir/grat_"+str(i)+"_"+str(j)+".npy",np.array(Image.open("output_dir/tile_"+str(7621*i)+"_"+str(7621*j)+".tif")))


