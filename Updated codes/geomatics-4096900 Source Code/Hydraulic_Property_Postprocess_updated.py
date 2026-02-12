from __future__ import division
import os
import pandas as pd
from time import perf_counter
from GeoFlood_Filename_Finder import cfg_finder

def process_folder(folder_path, comid_file):
    hydroprop_basetable = os.path.join(folder_path, "hydroprop-basetable.csv")
    networkmap_file = os.path.join(folder_path, "networkMapping.csv")
    handpropotxt = os.path.join(folder_path, "hydroprop-fulltable.csv")

    if not os.path.exists(hydroprop_basetable) or not os.path.exists(networkmap_file):
        print(f"Skipping folder {folder_path}: missing required files")
        return

    df_result = pd.read_csv(hydroprop_basetable)
    df_network = pd.read_csv(networkmap_file)
    df_n = pd.read_csv(comid_file)

    df_network = pd.merge(df_network, df_n, on='COMID')
    df_result = pd.merge(df_result, df_network, left_on='CatchId', right_on='HYDROID')
    df_result = df_result.drop('HYDROID', axis=1).rename(columns=lambda x: x.strip(" "))

    df_result['TopWidth (m)'] = df_result['SurfaceArea (m2)'] / df_result['LENGTHKM'] / 1000
    df_result['WettedPerimeter (m)'] = df_result['BedArea (m2)'] / df_result['LENGTHKM'] / 1000
    df_result['WetArea (m2)'] = df_result['Volume (m3)'] / df_result['LENGTHKM'] / 1000
    df_result['HydraulicRadius (m)'] = df_result['WetArea (m2)'] / df_result['WettedPerimeter (m)']
    df_result['HydraulicRadius (m)'].fillna(0, inplace=True)
    df_result['Discharge (m3s-1)'] = df_result['WetArea (m2)'] * \
                                     pow(df_result['HydraulicRadius (m)'], 2.0/3) * \
                                     pow(df_result['SLOPE'], 0.5) / df_result['Roughness']
    df_result['FloodAreaRatio'] = df_result['SurfaceArea (m2)'] / df_result['AREASQKM'] / 1000000

    if df_result['Discharge (m3s-1)'].isna().sum() == len(df_result):
        print(f"Empty DataFrame in folder {folder_path}, check hydroprop basetable and COMID_Roughness.csv")

    df_result.to_csv(handpropotxt, index=False)
    print(f"Saved: {handpropotxt}")


def main():
    geofloodHomeDir, projectName, DEM_name, chunk_status, input_fn, output_fn, hr_status = cfg_finder()

    comid_file = os.path.join(geofloodHomeDir, input_fn,
                              "Hydraulics", projectName, "COMID_Roughness.csv")

    channel_segments_dir = os.path.join(geofloodHomeDir, output_fn, "GIS", projectName, "ChannelSegments")

    for folder_name in os.listdir(channel_segments_dir):
        folder_path = os.path.join(channel_segments_dir, folder_name)
        if not os.path.isdir(folder_path):
            continue
        process_folder(folder_path, comid_file)


if __name__ == '__main__':
    t0 = perf_counter()
    main()
    t1 = perf_counter()
    print("Total time taken to postprocess hydraulic properties:", t1 - t0, "seconds")
