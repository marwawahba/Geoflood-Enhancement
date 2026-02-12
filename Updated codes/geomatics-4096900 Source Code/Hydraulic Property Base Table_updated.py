import os
import subprocess

hand_tif = r"C:\Users\Ahmed\Desktop\GeoFlood-master\GeoOutputs\GIS\my_project\dem_hand.tif"
slp_tif = r"C:\Users\Ahmed\Desktop\GeoFlood-master\GeoOutputs\GIS\my_project\dem_slp.tif"
stage_txt = r"C:\Users\Ahmed\Desktop\GeoFlood-master\GeoInputs\Hydraulics\my_project\stage.txt"

channel_segments_dir = r"C:\Users\Ahmed\Desktop\GeoFlood-master\GeoOutputs\GIS\my_project\ChannelSegments"
tauDEM_exe = r"C:\Users\Ahmed\Desktop\GeoFlood-master\TauDEM\CatchHydroGeo" 

for folder_name in os.listdir(channel_segments_dir):
    folder_path = os.path.join(channel_segments_dir, folder_name)
    if not os.path.isdir(folder_path):
        continue

    catch_tif = os.path.join(folder_path, "dem_segmentCatchment.tif")
    river_attr_txt = os.path.join(folder_path, "River_Attribute.txt")
    hydroprop_csv = os.path.join(folder_path, "hydroprop-basetable.csv")

    if not os.path.exists(catch_tif) or not os.path.exists(river_attr_txt):
        print(f"Skipping folder {folder_name}: missing required files")
        continue

    cmd = [
        "mpiexec", "-n", "4", tauDEM_exe,
        "-hand", hand_tif,
        "-catch", catch_tif,
        "-catchlist", river_attr_txt,
        "-slp", slp_tif,
        "-h", stage_txt,
        "-table", hydroprop_csv
    ]

    print(f"Processing folder: {folder_name}")
    try:
        subprocess.run(cmd, check=True)
        print(f"Saved hydroprop-basetable.csv in folder: {folder_name}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing folder {folder_name}: {e}")
