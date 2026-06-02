import xarray as xr
import os
import shutil
from pathlib import Path

SIMULATION_BASE = Path('/home/vicentev/scratch/vicentev/Simulation_V4')
INPUT_DIR = SIMULATION_BASE / 'RESTARTS' / 'R_2011'
OUTPUT_RESTART_DIR = SIMULATION_BASE / 'RESTARTS' / 'R_2012'

FILE_PREFIX = 'Tau_0_001_Ads_0_01_MC_0_2_Vel_Hx1_2_V4_2011' 

def process_simulations():
    OUTPUT_RESTART_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Starting batch process: {INPUT_DIR.name} -> NetCDF & {OUTPUT_RESTART_DIR.name} Restarts\n" + "="*60)

    for i in range(1, 13):
        zarr_name = f"{FILE_PREFIX}_P{i}.zarr"
        nc_name = f"{FILE_PREFIX}_P{i}.nc"
        restart_name = f"RESTART_2011_2012_P{i}.zarr"

        zarr_in_path = INPUT_DIR / zarr_name
        nc_out_path = INPUT_DIR / nc_name
        restart_out_path = OUTPUT_RESTART_DIR / restart_name
        
        if not zarr_in_path.exists():
            print(f"SKIPPING: P{i} (Source not found: {zarr_name})")
            continue

        print(f"\n--- Processing Part P{i} ---")
        
        try:
            with xr.open_zarr(zarr_in_path, chunks={}) as ds:
                
                print(f" 1/3 Writing NetCDF: {nc_name}...")
                ds.to_netcdf(nc_out_path, mode='w')
                
                print(f" 2/3 Extracting restart and saving to: {restart_name}...")
                restart = ds.isel(obs=slice(-1, None))
                restart.to_zarr(restart_out_path, mode='w')
            
            nc_success = nc_out_path.exists() and nc_out_path.stat().st_size > 0
            restart_success = restart_out_path.exists()

            if nc_success and restart_success:
                print(f" 3/3 VERIFIED: NetCDF ({nc_out_path.stat().st_size / 1e6:.2f} MB) and Restart Zarr created successfully.")
                print(f"     REMOVING original: {zarr_name}...")
                shutil.rmtree(zarr_in_path)
            else:
                print(f"     ERROR: Outputs missing or empty. Preserving original Zarr.")

        except Exception as e:
            print(f"     CRITICAL FAILURE on P{i}: {e}")
            print("     Source Zarr preserved. Moving to next file.")

    print("\n" + "="*60)
    print("All tasks processed.")
    print("="*60)

if __name__ == "__main__":
    process_simulations()