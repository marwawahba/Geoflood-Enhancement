from __future__ import division
from __future__ import print_function
import os
import numpy as np
from osgeo import gdal, ogr
from GeoFlood_Filename_Finder import cfg_finder
from time import perf_counter


def compute_segment_slope(m_array, z_array, step=10):
    """
    تقسيم المجرى كل 10 متر وحساب الميل بين اول واخر نقطة داخل كل جزء
    تجاهل الميل = صفر او سالب
    حساب الميديان من القيم الموجبة فقط
    """

    total_length = m_array[-1]  # الطول الحقيقي للمجرى داخل الـ SEGMENT
    positive_slopes = []        # قائمة الميل الموجب داخل أجزاء الـ 10 م

    # تقسيم كل 10 متر
    for start_m in np.arange(0, total_length, step):
        end_m = start_m + step

        # اختيار النقاط التي تقع داخل هذا الجزء
        idx = np.where((m_array >= start_m) & (m_array <= end_m))[0]

        if len(idx) < 2:
            continue

        i1 = idx[0]
        i2 = idx[-1]

        dz = z_array[i1] - z_array[i2]
        dx = m_array[i2] - m_array[i1]

        if dx <= 0:
            continue

        slope = dz / dx

        # تجاهل القيم السالبة أو الصفرية
        if slope > 0:
            positive_slopes.append(slope)

    if len(positive_slopes) == 0:
        return 0.00001  # لا يوجد ميل موجب — fallback

    # الميديان للقيم الموجبة فقط
    return np.median(positive_slopes)



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
    cat_layer.CreateField(ogr.FieldDefn("HYDROID", ogr.OFTInteger))
    dst_field = cat_layer.GetLayerDefn().GetFieldIndex("HYDROID")
    gdal.Polygonize(cat_rasterBand, mask, cat_layer, dst_field, [], callback=None)
    ds.Destroy()

    ds = driver.Open(segcat_shp, 1)
    cat_layer = ds.GetLayer()
    cat_layer.CreateField(ogr.FieldDefn('AreaSqKm', ogr.OFTReal))

    # حساب مساحة كل catchment
    for feat in cat_layer:
        geom = feat.GetGeometryRef()
        feat.SetField('AreaSqKm', float(geom.Area()) / 1000**2)
        cat_layer.SetFeature(feat)
        feat.Destroy()

    ac_iter = 0

    # ----------- حساب الميل الجديد داخل كل SEGMENT -------------
    for feature in layer:

        geom = feature.GetGeometryRef()
        feat_id = feature.GetField('HYDROID')
        length = feature.GetField('Length') / 1000

        points = geom.GetPoints()
        mx = np.array([p[0] for p in points])
        my = np.array([p[1] for p in points])

        px = ((mx - gt[0]) / gt[1]).astype(int)
        py = ((my - gt[3]) / gt[5]).astype(int)

        # حساب المسافة بين النقاط
        x_diff = np.diff(px) * gt[1]
        y_diff = np.diff(py) * gt[1]
        m_array = np.sqrt(x_diff**2 + y_diff**2)
        m_array = np.insert(m_array, 0, 0)
        m_array = np.cumsum(m_array)

        # ارتفاع النقاط
        z_array = rasterBand.ReadAsArray()[py, px].flatten()

        if np.sum(m_array) == 0:
            print("Empty Geometry Encountered")
            continue

        # ---------------- أهم نقطة: حساب الميل الجديد ----------------
        slope = compute_segment_slope(m_array, z_array, step=10)

        # ---------------- حساب المساحة ----------------
        cat_layer.SetAttributeFilter("HYDROID = {}".format(feat_id))
        area = 0
        for feat in cat_layer:
            area += feat.GetField("AreaSqKm")
            feat.Destroy()

        # تخزين القيم
        rafile.write("{} {} {} {}\n".format(feat_id, slope, length, area))

        ac_iter += 1
        print("Segment:", ac_iter)

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
            print(f"Skipping {folder_name}: shapefile missing")
            continue
        if not os.path.exists(segcatfn):
            print(f"Skipping {folder_name}: catchment raster missing")
            continue

        print("Processing:", folder_name)
        river_attribute_estimation(segment_shp, segcatfn, segcat_shp, demfn, attribute_txt)
        print("Saved:", folder_name)


if __name__ == '__main__':
    t0 = perf_counter()
    main()
    print("Time:", perf_counter() - t0)
