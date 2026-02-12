from __future__ import division
import os
import pandas as pd
from osgeo import ogr
from time import perf_counter 
from GeoFlood_Filename_Finder import cfg_finder

def network_mapping(cat_shp, seg_shp, map_csv):
    if not os.path.exists(cat_shp): 
        print("Catchment shapefile not found:", cat_shp)
        return
    if not os.path.exists(seg_shp):
        print("Segment shapefile not found:", seg_shp)
        return

    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(seg_shp, 0)
    layer = dataSource.GetLayer()
    hydroid_list = []
    comid_list = []

    for feature in layer:
        geom = feature.GetGeometryRef()
        pnt = geom.Centroid()
        cat_Source = driver.Open(cat_shp, 0)
        cat_layer = cat_Source.GetLayer()
        for cat_feat in cat_layer:
            if pnt.Intersects(cat_feat.GetGeometryRef()):
                hydroid = feature.GetField("HYDROID")
                comid = cat_feat.GetField("FEATUREID")
                hydroid_list.append(hydroid)
                comid_list.append(comid)
        cat_Source.Destroy()

    df = pd.DataFrame({"HYDROID": hydroid_list, "COMID": comid_list})
    df.to_csv(map_csv, index=False, columns=['HYDROID', 'COMID'])
    print("Saved network mapping:", map_csv)

def main():
    geofloodHomeDir, projectName, DEM_name, chunk_status, input_fn, output_fn, hr_status = cfg_finder()
    
    # البحث عن Catchment.shp في مجلد GeoInputs
    cat_shp = os.path.join(geofloodHomeDir, input_fn, "GIS", projectName, "Catchment.shp")
    if not os.path.exists(cat_shp):
        print(f"لم يتم العثور على Catchment.shp في: {cat_shp}")
        return
    
    # مجلد القنوات
    seg_base_dir = os.path.join(geofloodHomeDir, output_fn, "GIS", projectName, "ChannelSegments")
    
    # المرور على كل فولدر في ChannelSegments
    for folder_name in os.listdir(seg_base_dir):
        folder_path = os.path.join(seg_base_dir, folder_name)
        if not os.path.isdir(folder_path):
            continue

        seg_shp = os.path.join(folder_path, "dem_channelSegment.shp")
        map_csv = os.path.join(folder_path, "networkMapping.csv")
        
        print("Processing folder:", folder_name)
        network_mapping(cat_shp, seg_shp, map_csv)

if __name__ == '__main__':
    t0 = perf_counter()
    main()
    t1 = perf_counter()
    print("Total time taken:", t1-t0, "seconds")
