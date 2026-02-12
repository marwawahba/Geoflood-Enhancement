import os
import subprocess

# Paths
channel_segments_dir = r"C:\Users\Ahmed\Desktop\GeoFlood-master\GeoOutputs\GIS\my_project\ChannelSegments"
inunmap_base = r"C:\Users\Ahmed\Desktop\GeoFlood-master\GeoOutputs\Inundation\my_project"
hand_file = r"C:\Users\Ahmed\Desktop\GeoFlood-master\GeoOutputs\GIS\my_project\dem_hand.tif"

def runInunMap(hand_file, catch_file, forecast_nc, map_output):
    """
    Run InunMap to generate flood inundation map.
    """
    cmd = [
        "mpiexec", "-n", "4",
        r"C:\Users\Ahmed\Desktop\GeoFlood-master\TauDEM\InunMap",
        "-hand", hand_file,
        "-catch", catch_file,
        "-forecast", forecast_nc,
        "-mapfile", map_output
    ]
    print(f"\nRunning InunMap on {forecast_nc} ...")
    subprocess.run(cmd, check=True)
    print(f"Created inundation map: {map_output}")

# Forecast years
forecast_years = [5, 10, 25, 50, 100]

for folder_name in os.listdir(channel_segments_dir):
    folder_path = os.path.join(channel_segments_dir, folder_name)
    if not os.path.isdir(folder_path):
        continue

    # Create output folder
    inunmap_folder = os.path.join(inunmap_base, folder_name)
    os.makedirs(inunmap_folder, exist_ok=True)

    catch_file = os.path.join(folder_path, "dem_segmentCatchment.tif")

    for yr in forecast_years:
        forecast_nc = os.path.join(folder_path, f"forecast_{yr}yr.nc")
        map_output = os.path.join(inunmap_folder, f"dem_NWM_inunmap_{yr}yr.tif")

        if os.path.exists(catch_file) and os.path.exists(forecast_nc):
            runInunMap(hand_file, catch_file, forecast_nc, map_output)
        else:
            print(f"Missing file in folder {folder_name} for {yr}yr — skipped.")
