#This is a function that takes as arguments the resolution it is working with, and evolves a given array until the flow converges to a given level, or it runs out of steps.
import numpy as np
from matplotlib import pyplot as plt

def EvolveArray(path, resolution, convergence, steps, precipitation):
    
    res = resolution
    print("res = " + str(res))
    
    # convergence criterion parameter
    conv = convergence

    # Path to memory location for arrays of a particular resolution
    thisdir = path

    inputpath = thisdir + "res"+str(res)+"/"

    gratextents = [4,8]

    timestep = 0.05 #Much bigger than 0.1 and the integrator becomes unstable.

    # Keep precip constant to ensure that the tributaries get filled
    precip = precipitation

    number_of_steps=steps
    output_period = 10

    # Create ghost zone arrays. EW, NS, old and new. Contain only depth information.
    # Each graticule has a two element ghost zone. So the ghost zone for NS work must be 4*gratextents[0] x res*gratextents[1] and EW must be gratextents[0]*res x 4*gratextents[1]
    # It doesn't need to include the diagonals because they're irrelevant for the algo.
    # Note to self: This would have been much easier to implement as methods on a class. Live and learn. 
    # Second note to self: hooray for sunk cost fallacy
    ghostNSold = np.zeros([4*gratextents[0],res*gratextents[1]])
    ghostNSnew = np.zeros([4*gratextents[0],res*gratextents[1]])
    ghostEWold = np.zeros([res*gratextents[0],4*gratextents[1]])
    ghostEWnew = np.zeros([res*gratextents[0],4*gratextents[1]])

    # Allocate graticule memory
    # Mola, metric, depth, flowEW, flowNS, positive flow, negative flow, norm, delta
    graticule_space = np.zeros([res+4,res+4,9])

    # Fill in ghost zone depth numbers (initialization)
    # Imagine there are two domains, one next to the other, and not periodic boundaries (zero flow)
    # 1* 1 1 1 1 g1
    #          g2 2 2 2 2 2*
    # When domain 1 runs, it duplicates its first row to make 1*, then reads the first row of 2 to make g1. 
    # When it finishes, it writes its rightmost value to g2. So there's a "crossing of arms" between reading and writing. 
    # The only memory structure that will make any sense in the long term is one that is only ghost zones in the correct order.
    
    # For ghostNS, it only needs 6 rows, since the poles are redundant, but 8 rows simplifies addressing.
    # The 0th and 1st row are the interior 1 and -2th rows of domain latindex = 0, and so on .
    
    # Image four domains, next to each other, periodic boundaries
    # g1 1 1 1 1 g1     g3 3 3 3 3 g3     g1 1 1 1 1 g1
    #          g2 2 2 2 2 g2     g4 4 4 4 4 g4
    # Domain 1 reads 4 and 2 to get g1s, which requires some cyclic references. If the ghost zones are labeled in order, they are g4, g2, g1, g3, g2, g4, g3, g1, ...
    
    # In each cycle, it will read remote ghost zones, then write local ones. It's a beautiful thing. 
    
    # minalt allows the precipitation algorithm to normalize for heights. 

    minalt = 0
    totalalts = 0
    totalaltsmetric = 0
    for latindex in range(gratextents[0]):
        for lonindex in range(gratextents[1]):
            graticule_space = np.load(inputpath+"/array"+str(latindex)+str(lonindex)+".npy")
            ghostNSnew[4*latindex:4*latindex+2,res*lonindex:res*(lonindex+1)] = graticule_space[2:4,2:-2,2]
            ghostNSnew[4*latindex+2:4*latindex+4,res*lonindex:res*(lonindex+1)] = graticule_space[-4:-2,2:-2,2]
            ghostEWnew[res*latindex:res*(latindex+1),4*lonindex:4*lonindex+2] = graticule_space[2:-2,2:4,2]
            ghostEWnew[res*latindex:res*(latindex+1),4*lonindex+2:4*lonindex+4] = graticule_space[2:-2,-4:-2,2]
            
            # set minalt
            minalt = np.min([np.min(graticule_space[:,:,0]),minalt])
            totalalts += np.sum(graticule_space[2:-2,2:-2,0])
            totalaltsmetric += np.sum(graticule_space[2:-2,2:-2,0]*graticule_space[2:-2,2:-2,1])

    # rescale to set the minimum altitude to "0"
    totalalts += -res*res*gratextents[0]*gratextents[1]*minalt
    totalaltsmetric += -res*res*gratextents[0]*gratextents[1]*minalt*2/np.pi*1.0000136796523493

    # This is an array of numbers representing the total flow within the system. As a measure of convergence.
    total_flow = []
    total_rain = []
    total_water = []
    total_water_moment = []

    norm_factor = (res*res*gratextents[0]*gratextents[1])

    print("[step, depth, flow, evap, water alt]")
    #loop time steps
    for i in range(number_of_steps):
        if i%output_period == 1:
            print([i, total_water[-1]/norm_factor, total_flow[-1]/norm_factor, total_rain[-1]/norm_factor, total_water_moment[-1]/total_water[-1]])
        if i>2:
            if total_flow[-1]/total_flow[-3] > 1-conv:
                print("flow converged")
                #print("total flow:")
                #print(total_flow[-10])
                break
            if i==number_of_steps-1:
                print("flow didn't converge, try increasing number of steps")
                print("total flow:")
                print(total_flow[-10:])
                break
        

        # refresh ghost zones (begin the new time step with fresh space for old ghost zones)
        ghostEWold = ghostEWnew

        # refresh total evaporation to zero
        total_evap = 0

        # add a new term to total flow. At this point, doing EW and NS flows separately. 
        total_flow.append(0)
        total_water.append(0)
        total_rain.append(0)
        total_water_moment.append(0)

        # loop over graticules
        for latindex in range(gratextents[0]):
            for lonindex in range(gratextents[1]):
                # Load graticule
                graticule_space = np.load(inputpath+"/array"+str(latindex)+str(lonindex)+".npy")
            
                # update ghost zones from old ghost zone 
                # (which was the new ghost zone from the previous time step)
                graticule_space[2:-2,:2,2] = ghostEWold[res*latindex:res*(latindex+1),(4*lonindex-2)%(4*gratextents[1]):(4*lonindex-2)%(4*gratextents[1])+2]
                graticule_space[2:-2,-2:,2] = ghostEWold[res*latindex:res*(latindex+1),(4*lonindex+4)%(4*gratextents[1]):(4*lonindex+4)%(4*gratextents[1])+2]
            
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
                total_flow[-1] += np.sum(np.abs(graticule_space[2:-2,2:-2,8]))
                total_water[-1] += np.sum(graticule_space[2:-2,2:-2,2]*graticule_space[2:-2,2:-2,1])      
                total_water_moment[-1] += np.sum(graticule_space[2:-2,2:-2,0]*graticule_space[2:-2,2:-2,1]*graticule_space[2:-2,2:-2,2])
                
                # Note to self. Imagine an array like 
                # 0 0 1 0 0 0
                # Diff looks like
                # 0 1 -1 0 0 
                # Compute norm fraction for central 4
                # Diff+[:-1] - Diff-[1:]
                # 0 2 0 0 <--- norm = 0 0.5 0 0
                # Compute fate for central 4
                # (Diff+[1:] - Diff+[:-1] + Diff-[1:] - Diff-[:-1]) * norm
                # 1 0 0 0      0 -1 0 0     0 -1 0 0    0 0 1 0
                # 0.5 0 0 0    0 0.5 0 0    0 0.5 0 0   0 0 0.5 0
                
                # Determine correct evaporation level. Reusing layer 5 of graticule array.
                graticule_space[:,:,5] = np.minimum(graticule_space[:,:,2],precip)
                # Evaporate
                graticule_space[:,:,2] -= graticule_space[:,:,5]
                # Compute total volume by including metric, not including ghost zones
                total_evap += np.sum(graticule_space[2:-2,2:-2,5]*graticule_space[2:-2,2:-2,1])

                # update ghost zones (new)
                ghostEWnew[res*latindex:res*(latindex+1),4*lonindex:4*lonindex+2] = graticule_space[2:-2,2:4,2]
                ghostEWnew[res*latindex:res*(latindex+1),4*lonindex+2:4*lonindex+4] = graticule_space[2:-2,-4:-2,2]

                # save graticule
                np.save(inputpath+"/array"+str(latindex)+str(lonindex),graticule_space)



        # Do the same for NS
        # Cycle ghostzone. I hope this is assignation, not binding.
        ghostNSold = ghostNSnew
    
        total_flow.append(0)
        total_water.append(0)
        total_water_moment.append(0)
    
        # loop over graticules
        for latindex in range(gratextents[0]):
            for lonindex in range(gratextents[1]):
                # Load graticule
                graticule_space = np.load(inputpath+"/array"+str(latindex)+str(lonindex)+".npy")

                # update ghost zones from old ghost zone 
                # (which was the new ghost zone from the previous time step)
                if latindex == 0:
                    graticule_space[0,2:-2,2] = graticule_space[2,2:-2,2]
                    graticule_space[1,2:-2,2] = graticule_space[2,2:-2,2]
                else:
                    graticule_space[:2,2:-2,2] = ghostNSold[4*latindex-2:4*latindex,res*lonindex:res*(lonindex+1)]
                if latindex == gratextents[0]-1:
                    graticule_space[-1,2:-2,2] = graticule_space[-3,2:-2,2]
                    graticule_space[-2,2:-2,2] = graticule_space[-3,2:-2,2]
                else:
                    graticule_space[-2:,2:-2,2] = ghostNSold[4*latindex+4:4*latindex+6,res*lonindex:res*(lonindex+1)]
        
                # Rain including metric, including ghost zones
                # This is actually a volume
                rain_volume = total_evap*(graticule_space[:,:,0]-minalt)*graticule_space[:,:,1]/totalaltsmetric
            
                # volume of water
                total_rain[-1] += np.sum(rain_volume[2:-2,2:-2])
            
                #depth (convert by dividing by metric)
                graticule_space[:,:,2] += rain_volume/graticule_space[:,:,1]
            
                #volume
                total_water[-1] += np.sum(graticule_space[2:-2,2:-2,2]*graticule_space[2:-2,2:-2,1])
                total_water_moment[-1] += np.sum(graticule_space[2:-2,2:-2,2]*graticule_space[2:-2,2:-2,1]*graticule_space[2:-2,2:-2,0])

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

                # update ghost zones (new)
                ghostNSnew[4*latindex:4*latindex+2,res*lonindex:res*(lonindex+1)] = graticule_space[2:4,2:-2,2]
                ghostNSnew[4*latindex+2:4*latindex+4,res*lonindex:res*(lonindex+1)] = graticule_space[-4:-2,2:-2,2]

                # save graticule
                np.save(inputpath+"/array"+str(latindex)+str(lonindex),graticule_space)

        total_rain.append(total_evap)

    


