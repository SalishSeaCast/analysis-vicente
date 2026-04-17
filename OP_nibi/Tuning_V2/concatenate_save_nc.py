import xarray as xr
import os

# --- Configuration ---
path = '/home/vicentev/scratch/vicentev/Tuning_outputs/'

files_n = [
    'PBDEs_0112007_run_365_days_MC_0_2_Tau_0_0025_Ads_0_01_Vel_N.zarr',
    'PBDEs_0112007_run_365_days_MC_0_2_Tau_0_0025_Ads_0_01_Vel_Hx1_2.zarr',
    'PBDEs_0112007_run_365_days_MC_0_2_Tau_0_0025_Ads_0_1_Vel_N.zarr',
    'PBDEs_0112007_run_365_days_MC_0_2_Tau_0_0025_Ads_0_1_Vel_Hx1_2.zarr'
]

files_R = [
    'MC_0_2_Tau_0_0025_Ads_0_01_Vel_N_restarted.zarr',
    'MC_0_2_Tau_0_0025_Ads_0_01_Vel_Hx1_2_restarted.zarr',
    'MC_0_2_Tau_0_0025_Ads_0_1_Vel_N_restarted.zarr',
    'MC_0_2_Tau_0_0025_Ads_0_1_Vel_Hx1_2_restarted.zarr'
]

def concatenate_and_save():
    datasets = []
    for fn, fr in zip(files_n, files_R):
        # Load datasets
        ds1 = xr.open_zarr(os.path.join(path, fn))
        ds2 = xr.open_zarr(os.path.join(path, fr))
        
        # Concatenate
        ds_concat = xr.concat([ds1, ds2], dim='obs')  
        datasets.append(ds_concat)
        
        # Create the new filename with .nc extension
        out_filename = fn.replace('.zarr', '.nc')
        out_filepath = os.path.join(path, out_filename)
        
        # Save the concatenated dataset to NetCDF
        print(f"Saving {out_filename}...")
        ds_concat.to_netcdf(out_filepath)

    print("All files concatenated and saved successfully!")

if __name__ == "__main__":
    concatenate_and_save()