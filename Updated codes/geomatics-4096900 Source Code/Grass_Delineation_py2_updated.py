# -*- coding: utf-8 -*-
import os
import shutil
import sys
import subprocess
import grass.script as g
import grass.script.setup as gsetup

def segment_catchment_delineation(fdrfn, segshp, segcatfn):
    # Path to GRASS 7.6 batch file
    grass7bin = r'C:\Program Files\GRASS GIS 7.6\grass76.bat'

    # Get GISBASE
    startcmd = [grass7bin, '--config', 'path']
    p = subprocess.Popen(startcmd, shell=False,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        print("ERROR: Cannot find GRASS GIS 7 start script")
        sys.exit(-1)
    gisbase = out.strip().decode('utf-8')

    # GISDB directory
    gisdbdir = os.path.join(os.path.expanduser("~"), "Documents", "grassdata")
    locationGeonet = 'geonet'
    mapsetGeonet = 'geonetuser'
    grassGISlocation = os.path.join(gisdbdir, locationGeonet)

    # Remove existing location if exists
    if os.path.exists(grassGISlocation):
        shutil.rmtree(grassGISlocation)

    # Initialize GRASS
    gsetup.init(gisbase, gisdbdir, locationGeonet, 'PERMANENT')

    # Create location from raster (DEM/FDR)
    # If EPSG known, add epsg=<code>
    g.run_command('g.proj', georef=fdrfn, location=locationGeonet)

    # Create a new mapset
    gsetup.init(gisbase, gisdbdir, locationGeonet, 'PERMANENT')  # init PERMANENT first
    g.run_command('g.mapset', flags='c', mapset=mapsetGeonet,
                  location=locationGeonet, dbase=gisdbdir)
    gsetup.init(gisbase, gisdbdir, locationGeonet, mapsetGeonet)

    # Import raster and vector
    g.run_command('r.in.gdal', input=fdrfn, output='fdr', overwrite=True)
    g.run_command('g.region', raster='fdr')
    g.run_command('v.import', input=segshp, output='Segment')
    g.run_command('v.to.rast', input='Segment', use='attr',
                  output='stream', attribute_column='HYDROID')

    # Delineate basins
    g.run_command('r.stream.basins', overwrite=True,
                  direction='fdr', stream_rast='stream', basins='subbasins')

    # Export raster to GeoTIFF
    g.run_command('r.out.gdal', overwrite=True,
                  input='subbasins', type='Int16',
                  output=segcatfn, format='GTiff')

def main():
    # Paths
    geofloodHomeDir = r"C:\Users\Ahmed\Desktop\GeoFlood-master"
    DEM_name = "dem"
    projectName = "my_project"

    geofloodResultsDir = os.path.join(geofloodHomeDir, "GeoOutputs", "GIS", projectName)
    Name_path = os.path.join(geofloodResultsDir, DEM_name)
    fdrfn = Name_path + '_fdr.tif'

    # Folder with multiple subfolders containing shapefiles
    channel_segments_folder = os.path.join(geofloodResultsDir, "ChannelSegments")

    # Loop over each subfolder
    for subfolder in os.listdir(channel_segments_folder):
        folder_path = os.path.join(channel_segments_folder, subfolder)
        if os.path.isdir(folder_path):
            segshp = os.path.join(folder_path, "{}_channelSegment.shp".format(DEM_name))
            if not os.path.exists(segshp):
                print("Shapefile not found: {}, skipping...".format(segshp))
                continue
            segcatfn = os.path.join(folder_path, "{}_segmentCatchment.tif".format(DEM_name))
            print("Processing {} -> {}".format(segshp, segcatfn))
            segment_catchment_delineation(fdrfn, segshp, segcatfn)

    print("All segmentations completed successfully.")

if __name__ == '__main__':
    main()
