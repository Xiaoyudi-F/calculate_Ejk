import os
import re
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta
import time


# z500 data
data_root = r""  # data path
output_dir = r""

start_date = "20241201"
end_date = "20250310"

subfolders = ['00z']



def extract_z500(nc_file):
    """
    extract z500 from .nc files
    """
    try:
        with xr.open_dataset(nc_file) as ds:
            return ds["z500"].values
    except Exception as e:
        print(f"[ERROR] read NC failed {nc_file}: {e}")
        return None


def get_date_range(start_date, end_date):
    """
    Return a list of dates from start_date to end_date
    """
    start = datetime.strptime(start_date, "%Y%m%d")
    end = datetime.strptime(end_date, "%Y%m%d")
    return [start + timedelta(days=i) for i in range((end - start).days + 1)]


def extract_hour(filename):
    """
    Extract the forecast lead time from the filename
    """

    match = re.search(r"-(\d+)h-oper-fc_z500\.nc$", filename)
    return int(match.group(1)) if match else float('inf')


def process_grib_files(run_date, subfolder):
    """
    process files
    """
    date_str = run_date.strftime("%Y%m%d")
    folder_path = os.path.join(data_root, date_str, subfolder)
    if not os.path.exists(folder_path):
        print(f"[WARNING] The folder does not exist: {folder_path}")
        return None


    nc_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith("_z500.nc")]
    nc_files = sorted(nc_files, key=extract_hour)

    # 15-day forecast, 6-hour interval, should have 61 files
    if len(nc_files) != 61:
        print(f"[WARNING] {date_str}/{subfolder} Unexpected number of NC files: {len(nc_files)}")
        return None

    z500_list = []
    for nc_file in nc_files:
        data = extract_z500(nc_file)
        if data is None:
            return None
        z500_list.append(data)

    return np.array(z500_list)  # shape: (121, 721, 1440)


def get_forecast(valid_date, j, subfolder):
    run_date = valid_date - timedelta(days=j)
    # print("run_date:", run_date, "j:", j, "subfolder:", subfolder, sep="\n")
    print("run_date:", run_date, "j:", j)
    data = process_grib_files(run_date, subfolder)
    if data is None:
        return None
    index = j * 4
    if index >= data.shape[0]:
        print(f"[ERROR] 索引越界：run_date {run_date.strftime('%Y%m%d')} j={j} index={index}")
        return None
    return data[index, :, :]


def compute_weighted_rmse(diff, lats, region="Global"):
    """
    weighted root mean square error
    region: "Global", "NH" , "SH"
    """
    lats = np.linspace(90, -90, 721)

    if region == "NH":
        diff = diff[0:361, :]
        lats = lats[0:361]
    elif region == "SH":
        diff = diff[360:721, :]
        lats = lats[360:721]

    weights = np.cos(np.radians(lats))[:, np.newaxis]

    weighted_diff = (diff ** 2) * weights

    rmse_value = np.sqrt(np.sum(weighted_diff) / (np.sum(weights) * 1440))
    return rmse_value


def compute_rmse_across_days(valid_dates, subfolder, max_lead=15, region="Global"):
    max_lead = int(max_lead)
    ejk_sum = np.zeros((max_lead + 1, max_lead + 1))
    count = np.zeros((max_lead + 1, max_lead + 1))

    for valid_date in valid_dates:
        forecasts = {}
        valid_available = True
        for j in range(max_lead + 1):
            fc = get_forecast(valid_date, j, subfolder)
            if fc is None:
                valid_available = False
                print(f"[WARNING] valid_date {valid_date.strftime('%Y%m%d')} missing forecast for dya j={j} ")
                break
            forecasts[j] = fc
        if not valid_available:
            continue

        for j in range(max_lead):
            for k in range(j + 1, max_lead + 1):
                diff = forecasts[j] - forecasts[k]
                lat_dim = diff.shape[0]
                lats = np.linspace(90, -90, lat_dim)

                rmse = compute_weighted_rmse(diff, lats, region=region)
                ejk_sum[j, k] += rmse
                count[j, k] += 1

                print(f"[DEBUG] date: {valid_date.strftime('%Y%m%d')}, j={j}, k={k}, weighted RMSE={rmse}")

        print(f"[INFO]  {valid_date.strftime('%Y%m%d')} calculated.")

    with np.errstate(divide='ignore', invalid='ignore'):
        ejk_avg = np.where(count > 0, ejk_sum / count, np.nan)
    return ejk_avg


def main():
    os.makedirs(output_dir, exist_ok=True)
    valid_dates = get_date_range(start_date, end_date)

    # "Global", "NH", "SH"
    current_region = "Global"

    print(f"[INFO] Valid Date statistics: {len(valid_dates)}")


    print(f"[INFO] Start calculating the {current_region} weighted average Ejk (Lead=15) ...")
    ejk_avg = compute_rmse_across_days(valid_dates, subfolders[0], max_lead=15, region=current_region)

    output_filename = f" "
    output_path = os.path.join(output_dir, output_filename)
    pd.DataFrame(ejk_avg).to_csv(output_path, index=False, header=False)
    print(f"[INFO] Result saved: {output_path}")


if __name__ == "__main__":
    start = time.perf_counter()
    main()
    end = time.perf_counter()
    print("time consuming：", end - start)