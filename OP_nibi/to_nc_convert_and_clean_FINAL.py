import xarray as xr
import os
import shutil
from pathlib import Path

# --- Configuration ---
# Using Path objects for more robust cross-platform path handling
BASE_DIR = Path('/home/vicentev/scratch/vicentev/Simulation_V2/')

# Your specific list of files
ff = ['/home/vicentev/scratch/vicentev/Simulation_V2/PBDEs_0112007_run_365_days_Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_P1.zarr',
      '/home/vicentev/scratch/vicentev/Simulation_V2/PBDEs_0112007_run_365_days_Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_P2.zarr',
      '/home/vicentev/scratch/vicentev/Simulation_V2/PBDEs_0112007_run_365_days_Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_P3.zarr',
      '/home/vicentev/scratch/vicentev/Simulation_V2/PBDEs_0112007_run_365_days_Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_P4.zarr',
      '/home/vicentev/scratch/vicentev/Simulation_V2/PBDEs_0112007_run_365_days_Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_P5.zarr',
      '/home/vicentev/scratch/vicentev/Simulation_V2/PBDEs_0112007_run_365_days_Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_P6.zarr',
      '/home/vicentev/scratch/vicentev/Simulation_V2/PBDEs_0112007_run_365_days_Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_P7.zarr',
      '/home/vicentev/scratch/vicentev/Simulation_V2/PBDEs_0112007_run_365_days_Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_P8.zarr',
      '/home/vicentev/scratch/vicentev/Simulation_V2/PBDEs_0112007_run_365_days_Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_P9.zarr',
      '/home/vicentev/scratch/vicentev/Simulation_V2/PBDEs_0112007_run_365_days_Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_P10.zarr',
      '/home/vicentev/scratch/vicentev/Simulation_V2/PBDEs_0112007_run_365_days_Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_P11.zarr',
      '/home/vicentev/scratch/vicentev/Simulation_V2/PBDEs_0112007_run_365_days_Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_P12.zarr'
]

def safe_convert_and_delete():
    for zarr_name in ff:
        zarr_path = BASE_DIR / zarr_name
        nc_path = zarr_path.with_suffix('.nc')
        
        # 1. Check if source exists
        if not zarr_path.exists():
            print(f"SKIPPING: {zarr_name} (Source not found)")
            continue

        print(f"\n--- Starting: {zarr_name} ---")
        
        try:
            # 2. Perform the conversion
            # Using context manager 'with' ensures the file is closed after the block
            with xr.open_zarr(zarr_path, chunks={}) as ds:
                print(f"Writing NetCDF to: {nc_path.name}...")
                # mode='w' overwrites if an old failed attempt exists
                ds.to_netcdf(nc_path, mode='w')
            
            # 3. CRITICAL SAFETY CHECK
            # Check if .nc exists AND has a file size > 0 bytes
            if nc_path.exists() and nc_path.stat().st_size > 0:
                print(f"VERIFIED: {nc_path.name} created successfully ({nc_path.stat().st_size / 1e6:.2f} MB).")
                print(f"REMOVING: {zarr_name}...")
                shutil.rmtree(zarr_path)
            else:
                print(f"ERROR: {nc_path.name} is empty or missing. Preserving source Zarr.")

        except Exception as e:
            print(f"CRITICAL FAILURE on {zarr_name}: {e}")
            print("Source Zarr preserved. Moving to next file.")

    print("\n" + "="*30)
    print("All tasks processed.")
    print("="*30)

if __name__ == "__main__":
    safe_convert_and_delete()