# water-flow

Welcome to the Mars water flow simulator.
Casey Handmer 2018
License: CC-SA

System requirements:

\>120 gigs free HD space. Yes, 120 gigs.
\>8 gigs of ram.
\>At least one processor.


Instructions for use:

- Depends on numpy, colorsys, pandas, matplotlib, PIL, and whatever else yells when it runs. 
- Clone the git repo to your favored local directory.
- Download the raw MOLA data from https://drive.google.com/open?id=1swrkPXBYQd7yRPAP027BVSERCvDK6MK8 .
- Ensure that the folder RawMola is in your local directory. The total is about 6 gigs, it might take a while to download.
- Open Main.py and ensure that your local path is set, and that the global equivalent depth and precipitation parameter are set to something sensible. The default values seem to work.
- Run main with python 3. This will create subdirectories in which to save preliminary results, and then sequentially run the water flow simulation at higher and higher frequencies. On my c. 2011 laptop, it took about 3 days to run the whole thing. 
- When the data is run, use PlotArrays.py to generate and save pngs of each of 32 subsections of the planet, then assemble them/print them as you see fit. Experimentation on plotting the rich dataset is encouraged, and several demo color maps are given.


Future (unscheduled) planned improvements
- Increase resolution using HRSC dataset: https://astrogeology.usgs.gov/search/map/Mars/Topography/HRSC_MOLA_Blend/Mars_HRSC_MOLA_BlendDEM_Global_200mp_v2
- Improve precipitation model to include prevailing wind direction and snow.
- Provide a temperature sensitive vegetation model in a plot option. 
- Develop an erosion model.
Latest blog: https://caseyhandmer.wordpress.com/2018/11/29/mars-global-hydrology-at-full-mola-resolution/ 


What is this?

Mars Global Surveyor, a NASA mission flown to Mars, carried the Mars Orbital Laser Altimeter instrument. It produced a simply eye-watering dataset of altitudes all over Mars, with a resolution of about 460m. This simulation is a quick and dirty way to simulate what rivers and oceans might look like on Mars if it were terraformed.

This is the first decent-sized software project I've ever done in Python, so many of the design choices are pretty weird. The strangest one is the division of the planet up into 32 sections for sequential (not parallel) processing. If I had twenty lifetimes I might add parallel processing to increase speed, but the division was intended to allow the calculation to proceed within limited memory. The entire dataset is about 70 gigs at full resolution, which exceeds my humble laptop's capabilities.

The simulation models the flow of water downhill, collecting into rivers, lakes, and oceans. Much of Mars already has reasonably well-defined drainage basins, so much of it "makes sense". I also threw in a precipitation step to evaporate water and drop it back on the mountains, so that there would be a water cycle, and more importantly, plenty of water in tiny tributaries for me to plot. Because of aspects of the numerical method, convergence is not particularly physical. That is, the final dataset represents a situation where there is much more water in high altitude tributaries than lower rivers, which if played forward in time would result in a flood. It is worth remembering that the water cycle on Earth is never truly in equilibrium either!

When the code runs, it spits out a variety of numbers corresponding to depth, flow, evaporation, and average water altitude. Most of these should stay pretty steady over time, but flow does diminish as numerical artifacts of interpolation get "washed away". Average water altitude is a function of the degree of precipitation. I chose the precipitation level, of 1.5mm per time step, to give a hydrology with a familiar look - lots of rivers flowing to oceans. A higher rate of evaporation and precipitation simply short circuits the water cycle, increasing the average altitude of the water. This is something like mass snow and glaciation, with lots of water retained in the mountains. Conversely, a lower value for evaporation and precipitation leads to more water lurking at the lowest points in the landscape, the great northern ocean and the various impact basins.


Where is the island in the Hellas basin?

A popular trope of Mars literature, including Blue Mars (Kim Stanley Robinson) and Ilium (Dan Simmons) is that there is a rise in the middle of the Hellas basin that could be an island. Unfortunately this is an error on early (pre MOLA) maps because of a mass concentration the crust under Hellas. So there is no big island in the middle of Hellas. There could be some small ones, and there is a nice one in Gale Crater, which the Curiosity Rover is currently climbing up.
