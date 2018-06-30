# This script imports MOLA data from a local directory and creates the arrays at the appropriate resolution.

# These arrays are 2+5760/2^n (n = 1, 2, ..., 7) x 6ish, contain
# - MOLA topo data (double)
# - metric (sin latitude)
# - initial depth (initialize at constant GEL, or interpolate/decimate from other file)
# - final depth
# - flow EW
# - flow NS (strictly this is redundant but is a useful piece of data for plotting etc)


# This script also initializes an array containing ghost zone information, which is how adjacent sectors pass depth information back and forth.

# Path to memory location for arrays of a particular resolution

#from PIL import Image
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

path = "/home/handmer/Documents/Mars/water-flow/res45/"
