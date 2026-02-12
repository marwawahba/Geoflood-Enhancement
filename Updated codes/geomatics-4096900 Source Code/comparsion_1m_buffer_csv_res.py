import os
import rasterio
from rasterio.merge import merge
from rasterio import features
import numpy as np
import pandas as pd
from shapely.geometry import shape, mapping
from shapely.ops import unary_union
from scipy.ndimage import binary_fill_holes
import fiona
from fiona.crs import from_epsg
from openpyxl import Workbook

# ---------------------------------------
# Base settings (adjust paths if needed)
# ---------------------------------------
base_dir = r"C:\Users\Ahmed\Desktop\GeoFlood-master\GeoOutputs\Inundation\my_project"
ref_dir = os.path.join(base_dir, "Ras")
splits_summary_file = os.path.join(base_dir, "splits_summary.csv")

ref_files = ["5_yrs.tif", "10_yrs.tif", "25_yrs.tif", "50_yrs.tif", "100_yrs.tif"]
comp_files = [
    "dem_NWM_inunmap_5yr.tif",
    "dem_NWM_inunmap_10yr.tif",
    "dem_NWM_inunmap_25yr.tif",
    "dem_NWM_inunmap_50yr.tif",
    "dem_NWM_inunmap_100yr.tif"
]

splits_df = pd.read_csv(splits_summary_file)

out_merged_dir = os.path.join(base_dir, "merged_per_recurrence")
os.makedirs(out_merged_dir, exist_ok=True)
out_boundaries_dir = os.path.join(base_dir, "boundaries_per_recurrence")
os.makedirs(out_boundaries_dir, exist_ok=True)

# ---------------------------------------
# Helper: build mosaic for reference + computed rasters
# ---------------------------------------
def build_mosaic_for_recurrence(ref_path, comp_paths):
    datasets = []
    try:
        src_ref = rasterio.open(ref_path)
        datasets.append(src_ref)
    except Exception as e:
        for ds in datasets:
            ds.close()
        raise

    for p in comp_paths:
        try:
            ds = rasterio.open(p)
            datasets.append(ds)
        except Exception:
            pass

    if len(datasets) == 0:
        raise RuntimeError("No datasets to merge")

    mosaic, out_trans = merge(datasets)
    for ds in datasets:
        ds.close()

    mosaic = mosaic[0].astype('float32')
    meta = datasets[0].meta.copy() if len(datasets) > 0 else src_ref.meta.copy()
    meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[0],
        "width": mosaic.shape[1],
        "transform": out_trans,
        "dtype": "float32",
        "count": 1
    })
    nodata = src_ref.nodata if hasattr(src_ref, 'nodata') else None
    crs = src_ref.crs if hasattr(src_ref, 'crs') else None
    return mosaic, meta, out_trans, crs, nodata

# ---------------------------------------
# Build mosaic & boundary per recurrence
# ---------------------------------------
combined_boundaries = {}

for ref_file, comp_file in zip(ref_files, comp_files):
    print(f"\n=== Building mosaic & boundary for: {ref_file} ===")
    ref_path = os.path.join(ref_dir, ref_file)
    comp_paths = []
    for folder in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder)
        if not os.path.isdir(folder_path) or folder.lower() == "ras":
            continue
        candidate = os.path.join(folder_path, comp_file)
        if os.path.exists(candidate):
            comp_paths.append(candidate)

    if not os.path.exists(ref_path) and len(comp_paths) == 0:
        print(f"No reference or comp rasters found for {ref_file}, skipping.")
        continue
    if not os.path.exists(ref_path):
        print(f"Reference {ref_path} not found, using comps only.")

    try:
        mosaic, meta, transform, crs, nodata = build_mosaic_for_recurrence(ref_path, comp_paths)
    except Exception as e:
        print("Error building mosaic for", ref_file, ":", e)
        continue

    merged_out = os.path.join(out_merged_dir, f"combined_mosaic_{ref_file.replace('.tif','')}.tif")
    with rasterio.open(merged_out, "w", **meta) as dst:
        dst.write(mosaic, 1)
    print("Saved merged mosaic:", merged_out)

    mask = np.zeros_like(mosaic, dtype=np.uint8)
    mask[mosaic > 0] = 1
    mask_filled = binary_fill_holes(mask).astype(np.uint8)

    shapes_gen = features.shapes(mask_filled, mask=mask_filled, transform=transform)
    polygons = [shape(geom) for geom, val in shapes_gen]
    if len(polygons) == 0:
        print(f"No polygons extracted for {ref_file}, skipping boundary creation.")
        continue

    final_poly = unary_union(polygons).buffer(0)
    buffered_poly = final_poly.buffer(1)

    out_shp = os.path.join(out_boundaries_dir, f"combined_boundary_{ref_file.replace('.tif','')}_buffer5.shp")
    schema = {"geometry": "Polygon", "properties": {}}
    if crs is None:
        with fiona.open(out_shp, "w", driver="ESRI Shapefile", schema=schema) as dst:
            dst.write({"geometry": mapping(buffered_poly), "properties": {}})
    else:
        with fiona.open(out_shp, "w", driver="ESRI Shapefile", crs=from_epsg(crs.to_epsg()), schema=schema) as dst:
            dst.write({"geometry": mapping(buffered_poly), "properties": {}})

    print("Saved boundary shapefile:", out_shp)

    combined_boundaries[ref_file] = {
        "polygon": buffered_poly,
        "transform": transform,
        "meta": meta,
        "nodata": nodata
    }

# ---------------------------------------
# Compare per folder against boundary
# ---------------------------------------
all_results = []

for folder in os.listdir(base_dir):
    comp_dir = os.path.join(base_dir, folder)
    if not os.path.isdir(comp_dir) or folder.lower() == "ras":
        continue

    print(f"\n>>> Comparing folder: {folder}")
    results = []

    for ref_file, comp_file in zip(ref_files, comp_files):
        comp_path = os.path.join(comp_dir, comp_file)
        if not os.path.exists(comp_path):
            continue

        ref_path = os.path.join(ref_dir, ref_file)
        if not os.path.exists(ref_path):
            print(f"Reference {ref_path} missing, using boundary only if present.")

        with rasterio.open(ref_path) as src_ref, rasterio.open(comp_path) as src_comp:
            ref_arr = src_ref.read(1)
            comp_arr = src_comp.read(1)
            ref_nodata = src_ref.nodata
            comp_nodata = src_comp.nodata
            transform = src_ref.transform

        if ref_file not in combined_boundaries:
            print(f"Boundary for {ref_file} not found, skipping comparison.")
            continue

        poly = combined_boundaries[ref_file]["polygon"]
        mask_boundary = features.geometry_mask([poly], invert=True, transform=transform, out_shape=ref_arr.shape)

        ref_masked = np.where(mask_boundary, ref_arr, np.nan)
        comp_masked = np.where(mask_boundary, comp_arr, np.nan)

        if (ref_nodata is not None):
            ref_masked[(ref_masked == ref_nodata) & (~np.isnan(comp_masked))] = 0
        if (comp_nodata is not None):
            comp_masked[(comp_masked == comp_nodata) & (~np.isnan(ref_masked))] = 0

        valid_mask = ~((np.isnan(ref_masked)) & (np.isnan(comp_masked)))
        ref_final = ref_masked[valid_mask]
        comp_final = comp_masked[valid_mask]

        ref_bin = ref_final > 0
        comp_bin = comp_final > 0
        TP = np.logical_and(ref_bin, comp_bin).sum()
        TN = np.logical_and(~ref_bin, ~comp_bin).sum()
        FP = np.logical_and(~ref_bin, comp_bin).sum()
        FN = np.logical_and(ref_bin, ~comp_bin).sum()

        TPf, TNf, FPf, FNf = map(float, [TP, TN, FP, FN])
        total = TPf + TNf + FPf + FNf
        accuracy = (TPf + TNf) / total if total > 0 else np.nan
        precision = TPf / (TPf + FPf) if (TPf + FPf) > 0 else np.nan
        recall = TPf / (TPf + FNf) if (TPf + FNf) > 0 else np.nan
        f1 = (2 * precision * recall) / (precision + recall) if precision > 0 and recall > 0 else np.nan
        numerator = TPf * TNf - FPf * FNf
        denominator = np.sqrt((TPf + FPf)*(TPf + FNf)*(TNf + FPf)*(TNf + FNf))
        mcc = numerator / denominator if denominator != 0 else np.nan
        fm = np.sqrt(precision * recall) if precision > 0 and recall > 0 else np.nan
        csi = TPf / (TPf + FPf + FNf) if (TPf + FPf + FNf) > 0 else np.nan
        merror = (FPf + FNf) / total if total > 0 else np.nan
        bias = (TPf + FPf) / (TPf + FNf) if (TPf + FNf) > 0 else np.nan
        roc = recall
        RFN = FNf / (TPf + FNf) if (TPf + FNf) > 0 else np.nan
        RTN = TNf / (TNf + FPf) if (TNf + FPf) > 0 else np.nan
        REP = FPf / (TPf + FPf) if (TPf + FPf) > 0 else np.nan
        RTP = TPf / (TPf + FNf) if (TPf + FNf) > 0 else np.nan

        split_row = splits_df[splits_df['Folder'] == folder]
        if not split_row.empty:
            num_splits = split_row['Number_of_splits'].values[0]
            split_len = split_row['Split_length'].values[0]
        else:
            num_splits = np.nan
            split_len = np.nan

        results.append({
            "Ref_File": folder,
            "Comp_File": comp_file,
            "TP": TP, "TN": TN, "FP": FP, "FN": FN,
            "Accuracy": accuracy, "Precision": precision, "Recall": recall,
            "F1": f1, "MCC": mcc, "FM": fm, "CSI": csi, "MError": merror,
            "Bias": bias, "ROC": roc, "RFN": RFN, "RTN": RTN, "REP": REP, "RTP": RTP,
            "Number_of_splits": num_splits,
            "Split_length": split_len
        })

    output_csv = os.path.join(comp_dir, "comparison_metrics.csv")
    pd.DataFrame(results).to_csv(output_csv, index=False)
    print("Saved per-folder comparison:", output_csv)
    if results:
        all_results.append(pd.DataFrame(results))

# Merge results and create Excel pivot
if all_results:
    merged = pd.concat(all_results, ignore_index=True)
    merged_csv_cleaned = os.path.join(base_dir, "merged_comparison_metrics_cleaned.csv")
    merged.to_csv(merged_csv_cleaned, index=False)
    print("Saved merged CSV:", merged_csv_cleaned)

    metrics = ["TP", "TN", "FP", "FN", "Accuracy", "Precision", "Recall",
               "F1", "MCC", "FM", "CSI", "MError", "Bias", "ROC",
               "Number_of_splits", "Split_length"]

    output_excel = os.path.join(base_dir, "comparison_pivot_buffer.xlsx")
    with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
        for metric in metrics:
            if metric not in merged.columns:
                continue
            pivot_table = merged.pivot(index='Ref_File', columns='Comp_File', values=metric)
            try:
                pivot_table.to_excel(writer, sheet_name=metric)
            except Exception:
                sheet = metric[:31]
                pivot_table.to_excel(writer, sheet_name=sheet)
    print("Saved Excel pivot:", output_excel)
else:
    print("No comparison results were produced.")

# Merge splits into final Excel
splits_file = os.path.join(base_dir, "splits_summary.csv")
pivot_file = os.path.join(base_dir, "comparison_pivot_buffer.xlsx")
output_file = os.path.join(base_dir, "comparison_pivot_buffer_with_splits_sorted.xlsx")

splits_df = pd.read_csv(splits_file)
splits_cols = ['Number_of_splits', 'Split_length']

try:
    pivot_data = pd.read_excel(pivot_file, sheet_name=None, engine='openpyxl')
except Exception as e:
    print(f"Error reading Excel: {e}")
    pivot_data = {}

merged_sheets = {}
for sheet_name, df in pivot_data.items():
    if "Ref_File" not in df.columns:
        merged_sheets[sheet_name] = df
        continue

    df_merged = pd.merge(df, splits_df[['Folder'] + splits_cols], how='left',
                         left_on='Ref_File', right_on='Folder')
    if 'Folder' in df_merged.columns:
        df_merged.drop(columns=['Folder'], inplace=True)

    cols = df_merged.columns.tolist()
    new_order = ['Ref_File'] + splits_cols + [c for c in cols if c not in ['Ref_File'] + splits_cols]
    df_merged = df_merged[new_order]
    df_merged = df_merged.sort_values(by='Split_length', ascending=True)
    merged_sheets[sheet_name] = df_merged

with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
    for sheet_name, df in merged_sheets.items():
        df.to_excel(writer, index=False, sheet_name=sheet_name)

print(f"Final Excel created: {output_file}")
