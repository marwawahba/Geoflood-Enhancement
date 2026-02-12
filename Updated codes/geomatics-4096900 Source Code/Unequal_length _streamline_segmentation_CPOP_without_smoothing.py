import geopandas as gpd
import rasterio
import numpy as np
from shapely.geometry import LineString
import matplotlib.pyplot as plt
import os

from src import CPOP, LogCost, continuous_piecewise_linear_approximation, reconstruct_segmentation

def median_abs_deviation(arr):
    arr = np.asarray(arr)
    med = np.median(arr)
    return np.median(np.abs(arr - med))


flowline_path = r"C:\Users\Ahmed\Desktop\GeoFlood-master\GeoInputs\GIS\my_project\Flowline.shp"
dem_path      = r"C:\Users\Ahmed\Desktop\GeoFlood-master\GeoInputs\GIS\my_project\dem.tif"
output_folder = r"C:\Users\Ahmed\Desktop\GeoFlood-master\GeoOutputs\GIS\my_project\ChannelSegments\d_NO_smoothing_CPOP_realSlope"
os.makedirs(output_folder, exist_ok=True)

flow_gdf = gpd.read_file(flowline_path)
line = flow_gdf.geometry.iloc[0]
dem = rasterio.open(dem_path)

# ---------------------------
def sample_line(line, step=10):
    distances = np.arange(0, line.length + 1e-6, step)  # include end
    points = [line.interpolate(d) for d in distances]
    return points, distances

points, distances = sample_line(line, step=10)
coords = [(p.x, p.y) for p in points]
elev = np.array([val[0] for val in dem.sample(coords)])


y_vals = elev.copy()


line_length = line.length
if line_length < 500:
    beta = 2
    scale = 1
elif line_length < 2000:
    beta = 2
    scale = 1
else:
    beta = 10
    scale = 1

sigma = median_abs_deviation(y_vals)
changepoints = CPOP(y_vals, beta, sigma=sigma, verbose=True)

changepoints = np.concatenate(([0], changepoints, [len(y_vals)-1]))
phis = reconstruct_segmentation(y_vals, changepoints, sigma, beta, LogCost(scale))
approx = continuous_piecewise_linear_approximation(changepoints, phis)


plt.figure(figsize=(12,5))
plt.plot(distances, y_vals, label="Elevation")
plt.plot(distances, approx, label="Segment Approx")
plt.scatter(distances[changepoints], np.array(approx)[changepoints], color='red')
plt.title("Flowline Segmentation (No Smoothing)")
plt.xlabel("Distance along flowline (meters)")
plt.ylabel("Elevation (m)")
plt.legend()
plt.grid()
plt.savefig(os.path.join(output_folder, "elevation_segmentation_no_smoothing.png"),
            dpi=300, bbox_inches="tight")
plt.close()


segments = []
lengths = []
slope_real = []

for i in range(len(changepoints)-1):
    start = changepoints[i]
    end   = changepoints[i+1]
    seg_points = points[start:end+1]
    seg_line = LineString([(p.x, p.y) for p in seg_points])
    segments.append(seg_line)
    
    
    lengths.append(seg_line.length)

 
    start_pt = points[start]
    end_pt   = points[end]
    x1, y1 = start_pt.x, start_pt.y
    x2, y2 = end_pt.x, end_pt.y

    z1 = y_vals[start]
    z2 = y_vals[end]

    dz = z1-z2
    dx = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

    if dx == 0:
        slope_val = np.nan
    else:
        slope_val = dz / dx

    slope_real.append(slope_val)


seg_gdf = gpd.GeoDataFrame({
    "HYDROID": np.arange(1, len(segments)+1),
    "Length": lengths,
    "Slope_real": slope_real
}, geometry=segments, crs=flow_gdf.crs)


seg_output_shp = os.path.join(output_folder, "dem_channelSegment.shp")
seg_output_csv = os.path.join(output_folder, "dem_channelSegment.csv")
seg_gdf.to_file(seg_output_shp)
seg_gdf.drop(columns="geometry").to_csv(seg_output_csv, index=False)

print("Done! Segments saved at:", seg_output_shp)
print("CSV saved at:", seg_output_csv)
