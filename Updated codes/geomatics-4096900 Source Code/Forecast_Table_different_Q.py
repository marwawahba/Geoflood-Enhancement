
from __future__ import division
import sys, os, string, time, re, getopt, glob, shutil, math
import netCDF4
import numpy as np
import pandas as pd
import csv
import configparser
import argparse
import inspect
from datetime import datetime
from time import perf_counter
from GeoFlood_Filename_Finder import cfg_finder


# Read discharge (Q) columns from Excel and map to network
def readForecast(input_file, df_netmap):
    df_excel = pd.read_excel(input_file, engine='openpyxl')
    t = datetime.now().strftime('%Y%m%d_%H%M%S')
    init_t = t

    q_columns = [col for col in df_excel.columns if col.startswith('Q_')]
    if not q_columns:
        raise ValueError("No Q_ columns found in the Excel file!")

    dfs = {}
    for col in q_columns:
        df_temp = df_excel[['COMID', col]].rename(columns={col: 'Q'})
        df_temp['Recurrence'] = col.replace('Q_', '')
        df_nwm = pd.merge(df_netmap, df_temp, on='COMID')
        dfs[col] = {'timestamp': t, 'init_timestamp': init_t, 'df_nwm': df_nwm}
    return dfs


# Compute stage (H) via interpolation on hydroprop table
def forecastH(init_timestr, timestr, hp_input, stage_output, df_nwm, recurrence):
    hpdata = pd.read_csv(hp_input)

    comids = df_nwm['HYDROID'].values
    Qs = df_nwm['Q'].values

    h = np.zeros_like(Qs, dtype=float)

    for i in range(len(comids)):
        h_array = hpdata[hpdata['CatchId'] == comids[i]]['Stage'].values
        q_array = hpdata[hpdata['CatchId'] == comids[i]]['Discharge (m3s-1)'].values

        if len(h_array) > 0 and len(q_array) > 0:
            h[i] = np.interp(Qs[i], q_array, h_array, right=-9999)
        else:
            h[i] = -9999

    saveForecast(init_timestr, timestr, stage_output, df_nwm, h, recurrence)


# Save results to CSV and NetCDF
def saveForecast(init_timestr, timestr, stage_output, df_nwm, h, recurrence):
    csv_output = stage_output.replace('.nc', f'_{recurrence}.csv')
    nc_output = stage_output.replace('.nc', f'_{recurrence}.nc')

    df = pd.DataFrame({
        "COMID": df_nwm['COMID'].values,
        "HYDROID": df_nwm['HYDROID'].values,
        "Recurrence": recurrence,
        "Q": df_nwm['Q'].values,
        "H": h
    })

    df.to_csv(csv_output, index=False)
    print(f"CSV saved: {csv_output}")

    rootgrp = netCDF4.Dataset(nc_output, "w", format="NETCDF4")
    rootgrp.Subject = f'Inundation forecast for {recurrence}'
    rootgrp.Initialization_Timestamp = init_timestr
    rootgrp.Timestamp = timestr

    index = rootgrp.createDimension("index", len(df))
    comid_var = rootgrp.createVariable("COMID", "u4", ("index",))
    hydroid_var = rootgrp.createVariable("HYDROID", "u4", ("index",))
    q_var = rootgrp.createVariable("Q", "f4", ("index",))
    h_var = rootgrp.createVariable("H", "f4", ("index",))

    comid_var[:] = df['COMID'].values
    hydroid_var[:] = df['HYDROID'].values
    q_var[:] = df['Q'].values
    h_var[:] = df['H'].values

    rootgrp.close()
    print(f"NetCDF saved: {nc_output}")


# Main
if __name__ == '__main__':
    input_file = r'C:\Users\user\Desktop\GeoFlood-master\GeoInputs\NWM\my_project\Q_values.xlsx'

    geofloodHomeDir, projectName, DEM_name, chunk_status, input_fn, output_fn, hr_status = cfg_finder()

    base_dir = os.path.join(geofloodHomeDir, output_fn, "GIS", projectName, "ChannelSegments")

    print("\nProcessing all folders in ChannelSegments...\n")

    for folder in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder)
        if not os.path.isdir(folder_path):
            continue

        print(f"\n==============================")
        print(f" Folder: {folder}")
        print(f"==============================")

        hp_input = os.path.join(folder_path, "hydroprop-fulltable.csv")
        netmap_file = os.path.join(folder_path, "networkMapping.csv")

        if not os.path.exists(hp_input) or not os.path.exists(netmap_file):
            print("Skipping folder: required files not found.")
            continue

        df_netmap = pd.read_csv(netmap_file)
        stage_output = os.path.join(folder_path, "forecast.nc")

        forecast_data = readForecast(input_file, df_netmap)

        for recurrence, tobj in forecast_data.items():
            timestr = tobj['timestamp']
            init_timestr = tobj['init_timestamp']
            df_nwm = tobj['df_nwm']
            rec = recurrence.replace('Q_', '')

            print(f" Running recurrence = {rec}")
            forecastH(init_timestr, timestr, hp_input, stage_output, df_nwm, rec)

    print("\nAll folders processed successfully!")
