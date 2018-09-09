import matplotlib.pyplot as plt
from planetaryimage import PDS3Image
import numpy as np

NN=0

image = PDS3Image.open("dtm_msl_gc0911.da4.pds")
altsarray2 = image.image #altitudes in meters
altsarray = altsarray2[::2**NN,::2**NN]

timestep = 0.05

precip = 0.000

init_depth = 30.0

inputpath = "/home/handmer/Documents/Mars/water-flow/"

number_of_steps = 50000
output_period = 50

# Mola, metric, depth, flowEW, flowNS, positive flow, negative flow, norm, delta
graticule_space = np.zeros([altsarray.shape[0]+4,altsarray.shape[1]+4,9])

# Set up altitudes. Flat boundary to enforce zero flow out of the box.
graticule_space[2:-2,2:-2,0] = altsarray
graticule_space[2:-2,1,0] = graticule_space[2:-2,2,0]
graticule_space[2:-2,0,0] = graticule_space[2:-2,2,0]
graticule_space[2:-2,-2,0] = graticule_space[2:-2,-3,0]
graticule_space[2:-2,-1,0] = graticule_space[2:-2,-3,0]
graticule_space[1,:,0] = graticule_space[2,:,0]
graticule_space[0,:,0] = graticule_space[2,:,0]
graticule_space[-2,:,0] = graticule_space[-3,:,0]
graticule_space[-1,:,0] = graticule_space[-3,:,0]

# Metric
graticule_space[:,:,1] += 1.0

# Depth. Initialize to GED of 60m?
graticule_space[:,:,2] += init_depth

# Import existing file
#graticule_space = np.load(inputpath+"/gale_"+str(init_depth)+"_"+str(precip)+".npy")
if NN<5:
    #graticule_space_sub = np.load(inputpath+"/gale_"+str(init_depth)+"_"+str(precip)+"_"+str(NN+1)+".npy")
    graticule_space_sub = np.load(inputpath+"/gale_"+str(init_depth)+"_"+str(0.0)+"_"+str(NN+1)+".npy")
    c1 = np.array(graticule_space_sub[2:-2,2:-2,2].shape)-np.array(graticule_space[2:-2:2,2:-2:2,2].shape)
    graticule_space[2:-2:2,2:-2:2,2]=graticule_space_sub[2:-2-c1[0],2:-2-c1[1],2]
    c2 = np.array(graticule_space_sub[2:-2,2:-2,2].shape)-np.array(graticule_space[3:-2:2,2:-2:2,2].shape)
    graticule_space[3:-2:2,2:-2:2,2]=graticule_space_sub[2:-2-c2[0],2:-2-c2[1],2]
    c3 = np.array(graticule_space_sub[2:-2,2:-2,2].shape)-np.array(graticule_space[2:-2:2,3:-2:2,2].shape)
    graticule_space[2:-2:2,3:-2:2,2]=graticule_space_sub[2:-2-c3[0],2:-2-c3[1],2]
    c4 = np.array(graticule_space_sub[2:-2,2:-2,2].shape)-np.array(graticule_space[3:-2:2,3:-2:2,2].shape)
    graticule_space[3:-2:2,3:-2:2,2]=graticule_space_sub[2:-2-c4[0],2:-2-c4[1],2]

# Stuff for precipitation
minalt = np.min(graticule_space[:,:,0])
totalalts = np.sum(graticule_space[:,:,0]-minalt)
totalaltsmetric = np.sum((graticule_space[2:-2,2:-2,0]-minalt)*graticule_space[2:-2,2:-2,1])

# Array of numbers representing flow statistics over time
total_flow = []
total_rain = []
total_water = []
total_water_moment = []

norm_factor = altsarray.shape[0]*altsarray.shape[1]

total_evap = 0
rain_volume = 0

print("[step, depth, flow, evap, water alt]")
# loop time steps
for i in range(number_of_steps):
    # Tweak boundary to prevent flow across it
    graticule_space[2:-2,1,2] = graticule_space[2:-2,2,2]
    graticule_space[2:-2,0,2] = graticule_space[2:-2,2,2]
    graticule_space[2:-2,-2,2] = graticule_space[2:-2,-3,2]
    graticule_space[2:-2,-1,2] = graticule_space[2:-2,-3,2]
    graticule_space[1,:,2] = graticule_space[2,:,2]
    graticule_space[0,:,2] = graticule_space[2,:,2]
    graticule_space[-2,:,2] = graticule_space[-3,:,2]
    graticule_space[-1,:,2] = graticule_space[-3,:,2]
    
    # Compute flow for EW in place
    graticule_space[:,:-1,3] = timestep*np.diff(graticule_space[:,:,0] + graticule_space[:,:,2], axis = 1)

    # Update graticule depths
    # This must be normalized by the existing depth, or depth will go negative. 
    # Important to preserve conservatism.
    # Compute positive flow (right to left)
    graticule_space[:,:-1,5] = (0.5+0.5*np.sign(graticule_space[:,:-1,3]))*graticule_space[:,:-1,3] #positive, flows to the left, loss from right, add to left
    # Compute negative flow (left to right)
    graticule_space[:,:-1,6] = (0.5-0.5*np.sign(graticule_space[:,:-1,3]))*graticule_space[:,:-1,3] #negative, flows to the right, loss from left, add to right

    # Compute norm
    graticule_space[:,:,7] *= 0.0
    # add in volume term
    graticule_space[:,:-1,5] *= graticule_space[:,:-1,1]
    graticule_space[:,:-1,6] *= graticule_space[:,:-1,1]
    # Compute norm
    graticule_space[:,1:-1,7] = np.minimum(graticule_space[:,1:-1,2]*graticule_space[:,1:-1,1]/(10**-8 + graticule_space[:,:-2,5] - graticule_space[:,1:-1,6]),1.0)
    # Adjust flows for norm
    graticule_space[:,:-1,5] *= graticule_space[:,1:,7]
    graticule_space[:,:-1,6] *= graticule_space[:,:-1,7]

    # update flow, including a metric term that converts flows to volumes, scaled by latitude.
    graticule_space[2:-2,2:-2,8] = (graticule_space[2:-2,2:-2,5] - graticule_space[2:-2,1:-3,5] + graticule_space[2:-2,2:-2,6] - graticule_space[2:-2,1:-3,6])/graticule_space[2:-2,2:-2,1]

    # save delta
    graticule_space[:,:,3] *= 0.0
    graticule_space[2:-2,2:-2,3] = graticule_space[2:-2,2:-2,8]

    # update depth
    graticule_space[2:-2,2:-2,2] += graticule_space[2:-2,2:-2,8]

    # Save stats
    total_flow.append(0)
    total_water.append(0)
    total_water_moment.append(0)
    total_flow[-1] += np.sum(np.abs(graticule_space[2:-2,2:-2,8]))
    total_water[-1] += np.sum(graticule_space[2:-2,2:-2,2]*graticule_space[2:-2,2:-2,1])      
    total_water_moment[-1] += np.sum(graticule_space[2:-2,2:-2,0]*graticule_space[2:-2,2:-2,1]*graticule_space[2:-2,2:-2,2])

    # Determine correct evaporation level. Reusing layer 5 of graticule array.
    graticule_space[:,:,5] = np.minimum(graticule_space[:,:,2],precip)
    # Evaporate
    graticule_space[:,:,2] -= graticule_space[:,:,5]
    # Compute total volume by including metric, not including ghost zones
    total_evap = np.sum(graticule_space[2:-2,2:-2,5]*graticule_space[2:-2,2:-2,1])
    # Rain including metric, including ghost zones
    # This is actually a volume
    rain_volume = total_evap*(graticule_space[:,:,0]-minalt)*graticule_space[:,:,1]/totalaltsmetric
            
    # volume of water
    total_rain.append(0)
    total_rain[-1] += np.sum(rain_volume[2:-2,2:-2])
    
    #depth (convert by dividing by metric)
    graticule_space[:,:,2] += rain_volume/graticule_space[:,:,1]
            
    #volume
    #total_water[-1] += np.sum(graticule_space[2:-2,2:-2,2]*graticule_space[2:-2,2:-2,1])
    #total_water_moment[-1] += np.sum(graticule_space[2:-2,2:-2,2]*graticule_space[2:-2,2:-2,1]*graticule_space[2:-2,2:-2,0])

    
    # Do the same for NS flow
    # Tweak boundary to prevent flow across it
    graticule_space[2:-2,1,2] = graticule_space[2:-2,2,2]
    graticule_space[2:-2,0,2] = graticule_space[2:-2,2,2]
    graticule_space[2:-2,-2,2] = graticule_space[2:-2,-3,2]
    graticule_space[2:-2,-1,2] = graticule_space[2:-2,-3,2]
    graticule_space[1,:,2] = graticule_space[2,:,2]
    graticule_space[0,:,2] = graticule_space[2,:,2]
    graticule_space[-2,:,2] = graticule_space[-3,:,2]
    graticule_space[-1,:,2] = graticule_space[-3,:,2]

    # Compute flow for NS in place
    graticule_space[:-1,:,4] = timestep*np.diff(graticule_space[:,:,0] + graticule_space[:,:,2], axis = 0)

    # Update graticule depths
    # This must be normalized by the existing depth, or depth will go negative. 
    # Important to preserve conservatism.
    # Compute positive flow (bottom to top)
    graticule_space[:-1,:,5] = (0.5+0.5*np.sign(graticule_space[:-1,:,4]))*graticule_space[:-1,:,4] #positive, flows to the left, loss from right, add to left
    # Compute negative flow (top to bottom)
    graticule_space[:-1,:,6] = (0.5-0.5*np.sign(graticule_space[:-1,:,4]))*graticule_space[:-1,:,4] #negative, flows to the right, loss from left, add to right
    # Compute norm
    graticule_space[:,:,7] *= 0.0
    # Add in volume term
    graticule_space[:-1,:,5] *= graticule_space[1:,:,1]
    graticule_space[:-1,:,6] *= graticule_space[:-1,:,1]
    graticule_space[1:-1,:,7] = np.minimum(graticule_space[1:-1,:,2]*graticule_space[1:-1,:,1]/(10**-8 + graticule_space[:-2,:,5] - graticule_space[1:-1,:,6]),1.0)
    # Adjust flows for norm
    graticule_space[:-1,:,5] *= graticule_space[1:,:,7]
    graticule_space[:-1,:,6] *= graticule_space[:-1,:,7]
    
    # reset flow with normalized parts, including a metric term that converts flows to volumes, scaled by latitude.
    # update depth
    graticule_space[2:-2,2:-2,8] = (graticule_space[2:-2,2:-2,5] - graticule_space[1:-3,2:-2,5] + graticule_space[2:-2,2:-2,6] - graticule_space[1:-3,2:-2,6])/graticule_space[2:-2,2:-2,1]

    # save flow
    graticule_space[:,:,4] *= 0.0
    graticule_space[2:-2,2:-2,4] = graticule_space[2:-2,2:-2,8]
    
    graticule_space[2:-2,2:-2,2] += graticule_space[2:-2,2:-2,8]
    total_flow[-1] += np.sum(np.abs(graticule_space[2:-2,2:-2,8]))

    total_flow.append(0)
    total_water.append(0)
    total_water_moment.append(0)
    total_flow[-1] += np.sum(np.abs(graticule_space[2:-2,2:-2,8]))
    total_water[-1] += np.sum(graticule_space[2:-2,2:-2,2]*graticule_space[2:-2,2:-2,1])      
    total_water_moment[-1] += np.sum(graticule_space[2:-2,2:-2,0]*graticule_space[2:-2,2:-2,1]*graticule_space[2:-2,2:-2,2])

    if i%output_period == 0:
        print([i, total_water[-1]/norm_factor, total_flow[-1]/norm_factor, total_rain[-1]/norm_factor, total_water_moment[-1]/total_water[-1]])
        np.save(inputpath+"/gale_"+str(init_depth)+"_"+str(precip)+"_"+str(NN),graticule_space)
