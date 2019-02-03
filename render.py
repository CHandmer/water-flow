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
    latitude = (90.0-np.cumsum(0.*topo+180.0/(7*topo.shape[0]),axis=0)-(array_lat_index*topo.shape[0])/(7*topo.shape[0])*180.0)[:topo.shape[0]-4,:topo.shape[1]-4]
    altitude = (topo+depth)[2:-2,2:-2]
    angle = np.arctan2((topo+depth)[3:-1,2:-2]-(topo+depth)[1:-3,2:-2],2*21344.0/(14*topo.shape[0]))+latitude
    temperature = -25.0+45.0*np.cos(np.pi/180*angle)-altitude/200
    return temperature

# PNG writer helper function
def writeRGBArrayToPNG(filename,array):
    Image.fromarray((255*array).astype('uint8')).save(filename+".png")

def writeRGBArrayToJPG(filename,array):
    Image.fromarray((255*array).astype('uint8')).save(filename+".jpg")

# terrain renderer
# Depth data is plotted logarithmically.
# ocean 10**0-10**4
# greenery 10**-4 - 10**-2
# dry 10**-10 - 10**-9
def renderTerrain(topo,depth,vikingRGB):
    temp = getTemp(topo,depth,3,6)#-30.0
    logdepth = np.log10(depth[2:-2,2:-2])
    #gradientsigmoid = 0.5+0.5*np.tanh((np.diff(topo+depth,axis=0)[1:-2,2:-2]**2+np.diff(topo+depth,axis=1)[2:-2,1:-2]**2)**0.5-50.0)
    # There are five broad states: ice, rock, grass, forest, oceans
    # The transitions are more-or-less well defined.
    # These arrays give relative portions of each. 
    # Within each, other functions give actual color
    rock = 0.5-0.5*3*(logdepth+3.7)/np.sqrt(1.0+3.0**2*(logdepth+3.7)**2)
    ice = np.round((0.5-0.5*np.sign(temp+5)))*(1-rock) # hard transition below -5C, smoothing to rock
    ocean = np.round(0.5+0.5*np.sign(logdepth+1))*(1-ice) # hard transition above 10cm.
    forest = (0.5+0.5*8*(logdepth+3.02)/np.sqrt(1.0+8.0**2*(logdepth+3.02)**2))*(0.5+0.5*1*(temp+3.5)/np.sqrt(1.0+1.0**2*(temp+3.5)**2))*(1-ocean)
    grass = np.maximum(1.0-ocean-rock-ice-forest,0.0)

    # specify color functions, either arrays or scalars
    rock_color = vikingRGB
    ice_color = 0.9*np.array([1.0,1.0,1.0])
    ocean_depth = np.minimum(10.0**logdepth,10.0)#+0*8.*gradientsigmoid,10.0)
    ocean_color = cl.hsv_to_rgb(np.transpose([temp*0.0+0.56,0.25+0.16*ocean_depth-0.011*ocean_depth**2,0.97-0.022*ocean_depth-0.006*ocean_depth**2],axes=[1,2,0]))
    forest_depth = np.minimum(ocean_depth,0.01)
    forest_color = cl.hsv_to_rgb(np.transpose([temp*0.0+0.23,(30-temp)/700.+0.44+45.0*forest_depth,-(30-temp)/800.+temp*0.0+0.22],axes=[1,2,0]))
    grass_depth = np.minimum(ocean_depth,0.0003)
    grass_color = cl.hsv_to_rgb(np.transpose([-(30-temp)/400.+0.23,-(30-temp)/400.+0.48+624.4*grass_depth+433000.*grass_depth**2,temp*0.0+0.4],axes=[1,2,0]))

    output = rock[:,:,np.newaxis]*rock_color + ice[:,:,np.newaxis]*ice_color + ocean[:,:,np.newaxis]*ocean_color + forest[:,:,np.newaxis]*forest_color + grass[:,:,np.newaxis]*grass_color
    # add randomness
    output += 0.05*np.random.rand(output.shape[0],output.shape[1],output.shape[2])-0.025
    return output
 
def reliefTransform(grad):
    return np.sign(grad)*np.abs(grad)**0.5

def addRelief(topodepth,RGBarray,angle):
    # Implement add relief function to alter brightness.
    gradNS = reliefTransform(np.diff(topodepth,axis=0)[1:-2,2:-2]/(21344.0/(14*topodepth.shape[0])))
    gradEW = reliefTransform(np.diff(topodepth,axis=1)[2:-2,1:-2]/(21344.0/(14*topodepth.shape[0])))
    return np.minimum(np.maximum((1.0+(gradNS*np.cos(angle)+gradEW*np.sin(angle))/30.0)[:,:,np.newaxis]*RGBarray,0.0),1.0)

# Specify the subarray
lati = 3
loni = 6

# Get a test batch of data
data = np.load("array_"+str(lati)+"_"+str(loni)+".npy")
topo = data[:,:,0]
depth = data[:,:,2]

# Produce the full resolution version of the Viking data
if(False):
    viking = Image.open("Mars_Viking_Color.jpg")
    vikingdata = np.asarray(viking)
    viking.close()
    writeRGBArrayToPNG("Viking_"+str(lati)+"_"+str(loni),getVikingRGBvector(lati,loni,vikingdata,topo.shape))
    del vikingdata

colordata=renderTerrain(topo, depth, np.asarray(Image.open("Viking_"+str(lati)+"_"+str(loni)+".png"))/255.0)
writeRGBArrayToPNG("render_"+str(lati)+"_"+str(loni),colordata)
writeRGBArrayToPNG("relief_"+str(lati)+"_"+str(loni),addRelief(topo,colordata,0.25*np.pi))

if(False):
    writeRGBArrayToJPG("render_"+str(lati)+"_"+str(loni),colordata)
    writeRGBArrayToJPG("relief_"+str(lati)+"_"+str(loni),addRelief(topo,colordata,0.25*np.pi))
