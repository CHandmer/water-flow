import os

pix = 7621
deg = 25.714

outputdir = "gdal_dir"
outputname = "terraformed_7621_aa_drier"

os.system("mkdir "+outputdir)
os.system("mkdir "+outputdir+"2")

# Run gdal_translate to specify edges of each of the pre-rendered tiles
for i in range(14):
    for j in range(7):
        # Reference locs -180 90 (TL), 180 90 (TR), 180 -90 (BR)
        # -180+i*deg, 90-j*deg
        # -180+(i+1)*deg, 90-j*deg
        # -180+(i+1)*deg, 90-(j+1)*deg
        os.system("gdal_translate -of VRT -a_srs EPSG:4326 -gcp 0 0 "+str(-180+i*deg)+" "+str(90-j*deg)+" -gcp "+str(pix)+" 0 "+str(-180+(i+1)*deg)+" "+str(90-j*deg)+" -gcp "+str(pix)+" "+str(pix)+" "+str(-180+(i+1)*deg)+" "+str(90-(j+1)*deg)+" renders/relief_"+str(j)+"_"+str(i)+"_aa.png gdal_dir/relief_"+str(j)+"_"+str(i)+"_aa.vrt")

# Run gdalwarp to stretch each tile to the right shape
for i in range(14):
    for j in range(7):
        os.system("gdalwarp -of VRT -t_srs EPSG:4326 gdal_dir/relief_"+str(j)+"_"+str(i)+"_aa.vrt "+outputdir+"2/relief_"+str(j)+"_"+str(i)+"_aa_warp.vrt")

# gdalbuiltvrt combines all the .vrt files into one.
os.system("gdalbuildvrt "+outputname+".vrt "+outputdir+"2/*.vrt")

# gdal2tiles generates a directory + out.kml file for a superoverlay compatible with google earth
os.system("gdal2tiles.py -p geodetic -k -z 0-9 "+outputname+".vrt")

# Admire the size of your dir
os.system("du -hs "+outputname)
