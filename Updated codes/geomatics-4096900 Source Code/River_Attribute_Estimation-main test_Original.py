from __future__ import division
from __future__ import print_function
import os
import numpy as np
from scipy.stats import gmean
from osgeo import gdal, ogr
from GeoFlood_Filename_Finder import cfg_finder
from time import perf_counter

def river_attribute_estimation(segment_shp, segcatfn,
                               segcat_shp, burndemfn,
                               attribute_txt):
    if not os.path.exists(segment_shp):
        print("Shapefile not found:", segment_shp)
        return
    if not os.path.exists(segcatfn):
        print("Segment catchment raster not found:", segcatfn)
        return

    rafile = open(attribute_txt, "w")
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(segment_shp, 0)
    if dataSource is None:
        print("Failed to open shapefile:", segment_shp)
        return
    layer = dataSource.GetLayer()
    srs = layer.GetSpatialRef()
    featureCount = layer.GetFeatureCount()
    rafile.write(str(featureCount) + "\n")

    raster = gdal.Open(burndemfn)
    gt = raster.GetGeoTransform()
    rasterBand = raster.GetRasterBand(1)

    cat_raster = gdal.Open(segcatfn)
    cat_rasterBand = cat_raster.GetRasterBand(1)
    mask = cat_rasterBand.GetMaskBand()

    ds = driver.CreateDataSource(segcat_shp)
    cat_layer = ds.CreateLayer(segcat_shp, srs)
    field = ogr.FieldDefn("HYDROID", ogr.OFTInteger)
    cat_layer.CreateField(field)
    dst_field = cat_layer.GetLayerDefn().GetFieldIndex("HYDROID")
    gdal.Polygonize(cat_rasterBand, mask, cat_layer, dst_field, [], callback=None)
    ds.Destroy()

    ds = driver.Open(segcat_shp, 1)
    cat_layer = ds.GetLayer()
    cat_layer.CreateField(ogr.FieldDefn('AreaSqKm', ogr.OFTReal))

    for feat in cat_layer:
        geom = feat.GetGeometryRef()
        feat.SetField('AreaSqKm', float(geom.Area()) / 1000**2)
        cat_layer.SetFeature(feat)
        feat.Destroy()

    # Initialize arrays for slope correction
    z_array_du = []
    m_array_du = []
    slope_array_du = []
    mx_first = []
    mx_last = []
    ac_iter = 0
    print(f'Total Segments: {featureCount}')

    for feature in layer:
        geom = feature.GetGeometryRef()
        feat_id = feature.GetField('HYDROID')
        length = feature.GetField('Length') / 1000
        points = geom.GetPoints()
        mx = np.array([p[0] for p in points])
        my = np.array([p[1] for p in points])
        mx_first.append(mx[0])
        mx_last.append(mx[-1])
        px = ((mx - gt[0]) / gt[1]).astype(int)
        py = ((my - gt[3]) / gt[5]).astype(int)
        x_diff = np.diff(px) * gt[1]
        y_diff = np.diff(py) * gt[1]
        m_array = np.sqrt(x_diff**2 + y_diff**2)
        m_array = np.insert(m_array, 0, 0)
        m_array = np.cumsum(m_array)
        z_array = rasterBand.ReadAsArray()[py, px].flatten()

        if np.sum(m_array) == 0:
            print('Empty Geometry Encountered')
            continue

       
        L_10 = int(len(m_array) * 0.1)
        L_85 = int(len(m_array) * 0.85)
        slope = (z_array[L_10] - z_array[L_85]) / (m_array[L_85] - m_array[L_10])

      
        if (slope <= 0.0000001) and (ac_iter == 0):
            slope = 0.00001

        subtraction_iter = 1
        prev_check = 0
        gmean_check = 0
        
        if (slope <= 0.0) and (ac_iter >= 1):
            previous_reach_index = ac_iter - subtraction_iter
            while (slope <= 0.000001) and (mx_first[ac_iter] == mx_last[previous_reach_index]) and (subtraction_iter <= 3):
                previous_reach_index = ac_iter - subtraction_iter
                prior_length = []
                prior_array_length = []
                for i in range(subtraction_iter):
                    prior_length.append(m_array_du[previous_reach_index + i][-1])
                    prior_array_length.append(np.size(m_array_du[previous_reach_index + i]))

                z_array_test = np.concatenate(z_array_du[previous_reach_index:])
                m_array_test = np.concatenate(m_array_du[previous_reach_index:])
                
                for i in range(len(prior_length)):
                    if len(prior_length) < 2:
                        m_array_test[prior_array_length[i]:] = m_array_test[prior_array_length[i]:] + prior_length[0]
                    else:
                        if (i + 1) < len(prior_length):
                            m_array_test[prior_array_length[i + 1]:(prior_array_length[i] + prior_array_length[i + 1])] = \
                            m_array_test[prior_array_length[i + 1]:(prior_array_length[i] + prior_array_length[i + 1])] + prior_length[i + 1]
                        else:
                            m_array_test[np.sum(prior_array_length):] = m_array_test[np.sum(prior_array_length):] + prior_length[0]
                
                L_negative_10 = int(np.size(m_array_test) * 0.1)
                L_negative_85 = int(np.size(m_array_test) * 0.85)
                slope = (z_array_test[L_negative_10] - z_array_test[L_negative_85]) / (m_array_test[L_negative_85] - m_array_test[L_negative_10])
                if slope > 0:
                    prev_check = 1
                subtraction_iter += 1
            
            if (slope <= 0.000001) and (subtraction_iter > 3):
                slope_array_du_np = np.asarray(slope_array_du)
                slope = gmean(slope_array_du_np[slope_array_du_np > 0])
                gmean_check = 1

    
        if slope > 0.0:
            slope_array_du.append(slope)
        else:
            slope_array_du_np = np.asarray(slope_array_du)
            slope = gmean(slope_array_du_np[slope_array_du_np > 0])
            slope_array_du.append(slope)

     
        m_array_du.append(m_array)
        z_array_du.append(z_array)

     
        cat_layer.SetAttributeFilter("HYDROID = {}".format(feat_id))
        area = 0
        for feat in cat_layer:
            area += feat.GetField("AreaSqKm")
            feat.Destroy()

        rafile.write("{} {} {} {}\n".format(feat_id, slope, length, area))
        ac_iter += 1
        print("Segment: {}".format(ac_iter))
        feature.Destroy()

    rafile.close()
    dataSource.Destroy()
    ds.Destroy()

def main():
    geofloodHomeDir, projectName, DEM_name, chunk_status, input_fn, output_fn, hr_status = cfg_finder()
    geofloodResultsDir = os.path.join(geofloodHomeDir, output_fn, "GIS", projectName)
    demfn = os.path.join(geofloodHomeDir, input_fn, "GIS", projectName, DEM_name + ".tif")

    channel_segments_dir = os.path.join(geofloodResultsDir, "ChannelSegments")

    for folder_name in os.listdir(channel_segments_dir):
        folder_path = os.path.join(channel_segments_dir, folder_name)
        if not os.path.isdir(folder_path):
            continue

        segment_shp = os.path.join(folder_path, "dem_channelSegment.shp")
        segcatfn = os.path.join(folder_path, "dem_segmentCatchment.tif")
        segcat_shp = os.path.join(folder_path, "segmentCatchment.shp")
        attribute_txt = os.path.join(folder_path, "River_Attribute.txt")

        if not os.path.exists(segment_shp):
            print(f"Skipping folder {folder_name}: shapefile not found")
            continue
        if not os.path.exists(segcatfn):
            print(f"Skipping folder {folder_name}: catchment raster not found")
            continue

        print(f"Processing river attributes for folder: {folder_name}")
        river_attribute_estimation(segment_shp, segcatfn, segcat_shp, demfn, attribute_txt)
        print(f"Saved attributes in folder: {folder_name}")

if __name__ == '__main__':
    t0 = perf_counter()
    main()
    t1 = perf_counter()
    print("Total time taken: {:.2f} seconds".format(t1 - t0))
