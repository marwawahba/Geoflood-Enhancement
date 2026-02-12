from __future__ import division
from __future__ import print_function

import os
import numpy as np
import pandas as pd
from scipy.stats import gmean
from osgeo import gdal, ogr
from GeoFlood_Filename_Finder import cfg_finder
from time import perf_counter


def read_slope_csv(csv_path):
    """
    Reads slope values from dem_channelSegment.csv
    Column 1: HYDROID
    Column 3: Slope
    """
    df = pd.read_csv(csv_path)
    slope_dict = dict(zip(df.iloc[:, 0], df.iloc[:, 2]))
    return slope_dict


def river_attribute_estimation(segment_shp, segcatfn,
                               segcat_shp, burndemfn,
                               attribute_txt, slope_dict):

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

    gdal.Polygonize(cat_rasterBand, mask, cat_layer,
                    dst_field, [], callback=None)
    ds.Destroy()

    ds = driver.Open(segcat_shp, 1)
    cat_layer = ds.GetLayer()
    cat_layer.CreateField(ogr.FieldDefn('AreaSqKm', ogr.OFTReal))

    for feat in cat_layer:
        geom = feat.GetGeometryRef()
        feat.SetField('AreaSqKm', float(geom.Area()) / 1000**2)
        cat_layer.SetFeature(feat)
        feat.Destroy()

    print(f'Total Segments: {featureCount}')

    ac_iter = 0

    for feature in layer:
        geom = feature.GetGeometryRef()
        feat_id = feature.GetField('HYDROID')
        length = feature.GetField('Length') / 1000

        if feat_id not in slope_dict:
            print(f"Slope not found for HYDROID {feat_id}")
            feature.Destroy()
            continue

        slope = slope_dict[feat_id]

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

    geofloodHomeDir, projectName, DEM_name, chunk_status, \
    input_fn, output_fn, hr_status = cfg_finder()

    geofloodResultsDir = os.path.join(
        geofloodHomeDir, output_fn, "GIS", projectName
    )

    demfn = os.path.join(
        geofloodHomeDir, input_fn, "GIS",
        projectName, DEM_name + ".tif"
    )

    channel_segments_dir = os.path.join(
        geofloodResultsDir, "ChannelSegments"
    )

    target_folders = [
        "d_NO_smoothing_CPOP_realSlope",
        "d_smoothed_CPOP_realSlope"
    ]

    for folder_name in target_folders:
        folder_path = os.path.join(channel_segments_dir, folder_name)

        if not os.path.isdir(folder_path):
            print(f"Folder not found: {folder_name}")
            continue

        segment_shp = os.path.join(folder_path, "dem_channelSegment.shp")
        segcatfn = os.path.join(folder_path, "dem_segmentCatchment.tif")
        segcat_shp = os.path.join(folder_path, "segmentCatchment.shp")
        attribute_txt = os.path.join(folder_path, "River_Attribute.txt")
        csv_path = os.path.join(folder_path, "dem_channelSegment.csv")

        if not os.path.exists(csv_path):
            print(f"CSV file not found in {folder_name}")
            continue

        slope_dict = read_slope_csv(csv_path)

        print(f"Processing river attributes for folder: {folder_name}")

        river_attribute_estimation(
            segment_shp,
            segcatfn,
            segcat_shp,
            demfn,
            attribute_txt,
            slope_dict
        )

        print(f"Saved attributes in folder: {folder_name}")


if __name__ == '__main__':
    t0 = perf_counter()
    main()
    t1 = perf_counter()
    print("Total time taken: {:.2f} seconds".format(t1 - t0))
