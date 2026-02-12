
# -*- coding: utf-8 -*-
import os
import csv
from osgeo import ogr
import pandas as pd
import numpy as np
import shutil

# ---------------------------------------------------
# Function to calculate split count and average length
# ---------------------------------------------------
def get_split_info(shapefile):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds = driver.Open(shapefile, 0)
    if ds is None:
        return 0, 0.0

    layer = ds.GetLayer()
    total_length = 0.0
    num_splits = 0

    for feature in layer:
        geom = feature.GetGeometryRef()
        if geom is not None:
            total_length += geom.Length()
            num_splits += 1

    ds.Destroy()
    return num_splits, (total_length / num_splits if num_splits else 0.0)


# ---------------------------
# MAIN SCRIPT
# ---------------------------
def main():

    base_folder = r"C:\Users\Ahmed\Desktop\GeoFlood-master\GeoOutputs\GIS\my_project\ChannelSegments"

    # Output CSV (final name unified here)
    output_csv = os.path.join(base_folder, "splits_summary.csv")

    # ==========================================================
    # 1) Create initial CSV with (Folder, Number_of_splits, Split_length)
    # ==========================================================
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Folder', 'Number_of_splits', 'Split_length'])

        for folder in sorted(os.listdir(base_folder)):
            folder_path = os.path.join(base_folder, folder)
            if not os.path.isdir(folder_path):
                continue

            shp_file = os.path.join(folder_path, "dem_channelSegment.shp")
            if not os.path.exists(shp_file):
                print(f"Shapefile missing in: {folder}")
                continue

            n, L = get_split_info(shp_file)
            writer.writerow([folder, n, round(L, 2)])
            print(f"Processed {folder}: {n} splits, length {L:.2f}")

    print("\nInitial CSV created.\n")

    # ==========================================================
    # 2) Clean, filter, and sort the CSV
    # ==========================================================
    df = pd.read_csv(output_csv)

    special = [
        "d_NO_smoothing_CPOP_realSlope",
        "d_smoothed_CPOP_realSlope"
    ]

    # Create Numeric_Value column
    def extract_value(x):
        if x in special:
            return x
        if isinstance(x, str) and x.startswith("d_") and x[2:].isdigit():
            return int(x[2:])
        try:
            return int(x.replace("d_", ""))
        except:
            return x

    df["Numeric_Value"] = df["Folder"].apply(extract_value)

    # Mark special rows
    df["is_special"] = df["Folder"].isin(special)

    # Create a numeric-safe sorting column
    df["Numeric_Safe"] = pd.to_numeric(df["Numeric_Value"], errors="coerce")
    df.loc[df["is_special"], "Numeric_Safe"] = np.inf

    # Separate normal and special rows
    df_normal = df[~df["is_special"]].copy()
    df_special = df[df["is_special"]].copy()

    # Sort normal rows by Numeric_Safe
    df_normal = df_normal.sort_values("Numeric_Safe")

    # Drop duplicates based on Split_length (keep smallest split distance)
    df_normal = df_normal.drop_duplicates(subset=["Split_length"], keep="first")

    # Merge
    df_cleaned = pd.concat([df_normal, df_special], ignore_index=True)

    # Final sorting
    df_cleaned = df_cleaned.sort_values(["is_special", "Numeric_Safe"])

    # Save final cleaned CSV
    df_cleaned.to_csv(output_csv, index=False, encoding='utf-8')
    print("Cleaned CSV saved.")

    # ==========================================================
    # 3) Delete folders NOT present in the cleaned CSV
    # ==========================================================
    print("\nChecking folders to remove...\n")

    valid_folders = df_cleaned["Folder"].tolist()

    for folder in os.listdir(base_folder):
        folder_path = os.path.join(base_folder, folder)

        if os.path.isdir(folder_path):
            if folder not in valid_folders:
                print(f"Deleting folder: {folder_path}")
                shutil.rmtree(folder_path)
            else:
                print(f"Keeping folder: {folder_path}")

    print("\nFolder cleanup complete.\n")

    # ==========================================================
    # 4) Copy CSV to Inundation folder
    # ==========================================================
    dest_dir = r"C:\Users\Ahmed\Desktop\GeoFlood-master\GeoOutputs\Inundation\my_project"
    dest_csv = os.path.join(dest_dir, "splits_summary.csv")

    os.makedirs(dest_dir, exist_ok=True)
    shutil.copy(output_csv, dest_csv)

    print(f"CSV copied to: {dest_csv}")

    # ==========================================================
    # 5) Create new folders in Inundation based on cleaned CSV
    # ==========================================================
    print("\nCreating folders in Inundation...\n")

    for folder in valid_folders:
        new_folder = os.path.join(dest_dir, folder)
        os.makedirs(new_folder, exist_ok=True)
        print(f"Created: {new_folder}")

    print("\nAll destination folders created successfully.\n")


# ---------------------------
if __name__ == "__main__":
    main()
