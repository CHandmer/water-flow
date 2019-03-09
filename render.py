#This script imports some array of topo and depth data, then renders and exports it as a png at full resolution.

# The core of the script is a function that converts topo and depth data to a color in HSV space. Input also includes the viking data.

from PIL import Image
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as cl

Image.MAX_IMAGE_PIXELS = 1000000000

# This function returns Viking color at a given point
def getVikingRGB(array_lat_index, array_lon_index, lat_index, lon_index,vikingdata,datadims):
    vikingdims = vikingdata.shape
    viking_lat=((array_lat_index)*(datadims[0]-4)+lat_index)/(7*(datadims[0]-4))*(vikingdims[0]-1)
    viking_lon=((array_lon_index)*(datadims[1]-4)+lon_index)/(14*(datadims[1]-4))*(vikingdims[1]-1)
    viking_lat_int=int(np.floor(viking_lat))
    viking_lon_int=int(np.floor(viking_lon))
    viking_lat_lo = viking_lat-viking_lat_int
    viking_lon_lo = viking_lon-viking_lon_int
    viking_rgb_grid = (vikingdata[viking_lat_int,viking_lon_int]/255.0,
                       vikingdata[viking_lat_int,(viking_lon_int+1)%vikingdims[1]]/255.0,
                       vikingdata[viking_lat_int+1,viking_lon_int]/255.0,
                       vikingdata[viking_lat_int+1,(viking_lon_int+1)%vikingdims[1]]/255.0)
    viking_rgb = (viking_rgb_grid[0]*(1.-viking_lat_lo)*(1.-viking_lon_lo) + viking_rgb_grid[1]*(1.-viking_lat_lo)*viking_lon_lo + viking_rgb_grid[2]*viking_lat_lo*(1.-viking_lon_lo) + viking_rgb_grid[3]*viking_lat_lo*viking_lon_lo)
    return viking_rgb

def getVikingRGBvector(array_lat_index, array_lon_index,vikingdata,datadims):
    return np.array([[getVikingRGB(array_lat_index, array_lon_index, lat_index, lon_index,vikingdata,datadims) for lon_index in range(0,datadims[0]-4)] for lat_index in range(0,datadims[0]-4)])

# Function that computes temperature of a given point
# Factors altitude, latitude, and local slope
def getTemp(topo, depth, array_lat_index, array_lon_index):
    unitres=topo.shape[0]-4
    latitude = (90.0-np.cumsum(0.*topo[2:-2,2:-2]+180.0/(7*unitres),axis=0)-(array_lat_index*unitres)/(7*unitres)*180.0)
    altitude = (topo+depth)[2:-2,2:-2]
    # Discount aspect angle by 0.5 since it is diluted by dawn and dusk.
    angle = 0.1*180.0/np.pi*np.arctan2((topo+depth)[3:-1,2:-2]-(topo+depth)[1:-3,2:-2],2*21344.0/(14*topo.shape[0]))+latitude
    # Decrease lapse rate to compensate for Mars' topographic variation. Increase peak average temp to 35C.
    #temperature = -25.0+60.0*np.cos(np.pi/180*angle)-altitude/400
    temperature = -25.0+60.0*np.cos(np.pi/180*angle)-altitude/200
    return temperature

# PNG writer helper function
def writeRGBArrayToPNG(filename,array):
    Image.fromarray((255*array).astype('uint8')).save(filename+".png")

def writeRGBArrayToJPG(filename,array):
    Image.fromarray((255*array).astype('uint8')).save(filename+".jpg")

def HSV2RGBhelper(outputarray, inputarray, mask):
    # Loop through by rows to save memory issue.
    for i in range(len(outputarray)):
        outputarray[i] += mask[i,:,np.newaxis]*cl.hsv_to_rgb(inputarray[i])
    #return cl.hsv_to_rgb(array)

# terrain renderer
# Depth data is plotted logarithmically.
# ocean 10**0-10**4
# greenery 10**-4 - 10**-2
# dry 10**-10 - 10**-9
def renderTerrain(topo,depth,vikingRGB,lati,loni):
    temp = getTemp(topo,depth,lati,loni)
    logdepth = np.log10(depth[2:-2,2:-2]+np.finfo(float).eps)
    #gradientsigmoid = 0.5+0.5*np.tanh((np.diff(topo+depth,axis=0)[1:-2,2:-2]**2+np.diff(topo+depth,axis=1)[2:-2,1:-2]**2)**0.5-50.0)
    # There are five broad states: ice, rock, grass, forest, oceans
    # The transitions are more-or-less well defined.
    # These arrays give relative portions of each. 
    # Within each, other functions give actual color

    # Original definitions
#    rock = 0.5-0.5*3*(logdepth+3.7)/np.sqrt(1.0+3.0**2*(logdepth+3.7)**2)
#    ice = np.round((0.5-0.5*np.sign(temp+5)))*(1-rock) # hard transition below -5C, smoothing to rock
#    ocean = np.round(0.5+0.5*np.sign(logdepth+1))*(1-ice) # hard transition above 10cm.
#    forest = (0.5+0.5*8*(logdepth+3.0)/np.sqrt(1.0+8.0**2*(logdepth+3.0)**2))*(0.5+0.5*1*(temp+3.5)/np.sqrt(1.0+1.0**2*(temp+3.5)**2))*(1-ocean)
#    grass = np.maximum(1.0-ocean-rock-ice-forest,0.0)

# New definitions
    rock = 0.5-0.5*3*(logdepth+3.7)/np.sqrt(1.0+3.0**2*(logdepth+3.7)**2)
    ice = np.round((0.5-0.5*np.sign(temp+5)))*(1-rock) # hard transition below -5C, smoothing to rock
    ocean = np.round(0.5+0.5*np.sign(logdepth+1))*(1-ice) # hard transition above 10cm.
    forest = (0.5+0.5*8*(logdepth+2.85)/np.sqrt(1.0+8.0**2*(logdepth+2.85)**2))*(0.5+0.5*1*(temp+1.5)/np.sqrt(1.0+1.0**2*(temp+1.5)**2))*(1-ocean)
    grass = np.maximum(1.0-ocean-rock-ice-forest,0.0)

    #Make tundra redder. Cut trees off in the hottest areas?
    
    # specify color functions, either arrays or scalars
    output = rock[:,:,np.newaxis]*vikingRGB
    del rock
    ice_color = 0.9*np.array([1.0,1.0,1.0])
    output += ice[:,:,np.newaxis]*ice_color
    del ice, ice_color
    ocean_depth = np.minimum(10.0**logdepth,10.0)#+0*8.*gradientsigmoid,10.0)
    del logdepth
    ocean_color_hsv = np.stack([temp*0.0+0.56,0.3+0.15*ocean_depth-0.011*ocean_depth**2,0.8-0.02*ocean_depth-0.004*ocean_depth**2],axis=-1)
    HSV2RGBhelper(output,ocean_color_hsv,ocean)
    del ocean, ocean_color_hsv
    forest_depth = np.minimum(ocean_depth,0.01)
    forest_color_hsv = np.stack([temp*0.0+0.23,(30-temp)/700.+0.44+45.0*forest_depth,-(30-temp)/800.+temp*0.0+0.22],axis=-1)
    HSV2RGBhelper(output,forest_color_hsv,forest)
    del forest, forest_color_hsv
    grass_depth = np.minimum(ocean_depth,0.0003)
    del ocean_depth
    #grass_color_hsv = np.stack([-(30-temp)/400.+0.23,-(30-temp)/400.+0.48+624.4*grass_depth+433000.*grass_depth**2,temp*0.0+0.4],axis=-1)
    grass_color_hsv = np.stack([-(30-temp)/300.+0.23,-(30-temp)/400.+0.48+624.4*grass_depth+433000.*grass_depth**2,temp*0.0+0.4],axis=-1)
    del grass_depth
    HSV2RGBhelper(output,grass_color_hsv,grass)
    del grass, grass_color_hsv

    #output = rock[:,:,np.newaxis]*rock_color + ice[:,:,np.newaxis]*ice_color + ocean[:,:,np.newaxis]*ocean_color + forest[:,:,np.newaxis]*forest_color + grass[:,:,np.newaxis]*grass_color
    # add randomness
    #output += 0.05*np.random.rand(output.shape[0],output.shape[1],output.shape[2])-0.025
    #print(np.sum(output))
    return output

def reliefTransform(grad):
    return np.sign(grad)*np.abs(grad)**0.5

def addRelief(topodepth,RGBarray,angle):
    # Implement add relief function to alter brightness.
    gradNS = reliefTransform(np.diff(topodepth,axis=0)[1:-2,2:-2]/(21344.0/(14*topodepth.shape[0])))
    gradEW = reliefTransform(np.diff(topodepth,axis=1)[2:-2,1:-2]/(21344.0/(14*topodepth.shape[0])))
    return np.minimum(np.maximum((1.0+(gradNS*np.cos(angle)+gradEW*np.sin(angle))/30.0)[:,:,np.newaxis]*RGBarray,0.0),1.0)

def antiAlias(image):
    imagesize=image.size
    return image.resize((2*imagesize[0],2*imagesize[1])).resize(imagesize,resample=Image.ANTIALIAS)

def renderImage(latindex,lonindex):
    
    # Specify the subarray
    lati = latindex
    loni = lonindex

    # Get a test batch of data
    print("importing data for lat: " + str(lati) + " and lon: " + str(loni))
    data = np.load("res7621_compressed/array_"+str(lati)+"_"+str(loni)+".npy")
    topo = data[:,:,0]
    depth = data[:,:,1]
    del data

    # Produce the full resolution version of the Viking data
    if(False):
        print("producing the viking data map")
        viking = Image.open("Mars_Viking_Color.jpg")
        vikingdata = np.asarray(viking)
        viking.close()
        writeRGBArrayToPNG("renders/Viking_"+str(lati)+"_"+str(loni),getVikingRGBvector(lati,loni,vikingdata,topo.shape))
        del vikingdata

    print("producing color data for render")
    colordata=renderTerrain2(topo, depth, np.asarray(Image.open("renders/Viking_"+str(lati)+"_"+str(loni)+".png"))/255.0,lati,loni)
    print("saving images")
    writeRGBArrayToPNG("renders/render_"+str(lati)+"_"+str(loni),colordata)
    writeRGBArrayToPNG("renders/relief_"+str(lati)+"_"+str(loni),addRelief(topo,colordata,0.25*np.pi))
    antiAlias(Image.open("renders/relief_"+str(lati)+"_"+str(loni)+".png")).save("renders/relief_"+str(lati)+"_"+str(loni)+"_aa.png")

    if(True):
        print("saving jpgs")
        writeRGBArrayToJPG("renders/render_"+str(lati)+"_"+str(loni),colordata)
        writeRGBArrayToJPG("renders/relief_"+str(lati)+"_"+str(loni),addRelief(topo,colordata,0.25*np.pi))
        antiAlias(Image.open("renders/relief_"+str(lati)+"_"+str(loni)+".png")).save("renders/relief_"+str(lati)+"_"+str(loni)+"_aa.jpg")


#for i in range(7):
#    for j in range(14):
#        renderImage(i,j)

renderImage(3,6)
#renderImage(2,1)

#This is a toy function that produces an image of the colormap as a function of temp and depth.
def colorMap(pix):
    topo = np.zeros(((pix+4),(pix+4)))+5000
    depth = np.array([[10.0**(-6.+10.*i/pix) for i in range(pix+4)] for j in range(pix+4)])
    rockcolor = np.zeros((pix,pix,3))+np.mean(np.mean(np.asarray(Image.open("Viking_3_6.png"))/255.0,axis=0),axis=0)

    print("producing color data for render")
    colordata=renderTerrain(topo, depth, rockcolor, 1, 0)
    print("saving images")
    writeRGBArrayToPNG("colorMap",colordata)

    if(True):
        print("saving jpgs")
        writeRGBArrayToJPG("colorMap",colordata)

