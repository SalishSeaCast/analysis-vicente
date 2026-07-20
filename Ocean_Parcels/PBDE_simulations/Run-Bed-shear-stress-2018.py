import xarray as xr
import numpy as np
import pandas as pd
import glob
import os
import dask
from dask.diagnostics import ProgressBar

def get_file_list(start, end, path, prefix_type):
    """Efficiently builds a list of NEMO output files using pandas date_range."""
    dates = pd.date_range(start, end, freq='D')
    file_list = []
    for date in dates:
        folder = date.strftime("%d%b%y").lower()
        file_str = f'SalishSea_1h_{date.strftime("%Y%m%d")}_{date.strftime("%Y%m%d")}_{prefix_type}*'
        search_path = os.path.join(path, folder, file_str)
        files = glob.glob(search_path)
        if files:
            file_list.extend(files)
    return sorted(file_list)

# Preprocess functions to strip out useless metadata/coordinates at load time
def prep_U(ds): return ds[['vozocrtx']]
def prep_V(ds): return ds[['vomecrty']]
def prep_T(ds): return ds[['e3t']]

def main():
    # --- 1. SETUP PATHS & FILES ---
    begin = pd.to_datetime("2018-01-01")
    end = pd.to_datetime("2018-12-31")
    
    path_days_h = '/results2/SalishSea/nowcast-green.202111/'
    filepaths_U = get_file_list(begin, end, path_days_h, 'grid_U.nc')
    filepaths_V = get_file_list(begin, end, path_days_h, 'grid_V.nc')
    filepaths_e3t = get_file_list(begin, end, path_days_h, 'grid_T.nc')

    # --- 2. LOAD MESH MASK & BATHYMETRY ---
    path_bat = '/ocean/vvalenzuela/MOAD/grid2/mesh_mask202108_TDV.nc'
    
    # Load just the mbathy array directly into memory to prevent Dask overhead
    with xr.open_dataset(path_bat) as bat_file:
        mbathy = bat_file['mbathy'].isel(t=0).load()

    # Highly optimized index clamping (clip avoids expensive .where copies)
    bottom_idx = (mbathy - 1).clip(min=0).astype(np.int16)
    
    # Create an in-memory boolean mask for land (True = ocean, False = land)
    ocean_mask = (mbathy > 0)

    # --- 3. LAZY LOAD ALL DATA (Aggressively Stripped) ---
    chunks = {'time_counter': 24}

    print("Building lean lazy datasets...")
    ds_U = xr.open_mfdataset(filepaths_U, chunks=chunks, parallel=True, preprocess=prep_U, coords='minimal', compat='override', engine='netcdf4')
    ds_V = xr.open_mfdataset(filepaths_V, chunks=chunks, parallel=True, preprocess=prep_V, coords='minimal', compat='override', engine='netcdf4')
    ds_e3t = xr.open_mfdataset(filepaths_e3t, chunks=chunks, parallel=True, preprocess=prep_T, coords='minimal', compat='override', engine='netcdf4')

    # --- 4. VECTORIZED COMPUTATION (Forced to float32) ---
    print("Setting up computational graph...")
    
    # Select bottom layer and explicitly cast to float32 to halve RAM usage
    vel_U = ds_U['vozocrtx'].isel(depthu=bottom_idx).astype(np.float32)
    vel_V = ds_V['vomecrty'].isel(depthv=bottom_idx).astype(np.float32)
    e3t = ds_e3t['e3t'].isel(deptht=bottom_idx).astype(np.float32)

    z = e3t / 2.0
    z_star = 0.07
    k = 0.42

    U_left = vel_U.shift(x=1)
    V_bottom = vel_V.shift(y=1)

    U_horizontal_2 = 0.25 * (vel_U + U_left)**2
    V_horizontal_2 = 0.25 * (vel_V + V_bottom)**2

    vel_horizontal = np.sqrt(U_horizontal_2 + V_horizontal_2)

    # Calculate and mask in one step
    u_star = ((vel_horizontal * k) / np.log(z / z_star)).where(ocean_mask)
    tau = ((u_star**2) * 1024.0).where(ocean_mask)

    # --- 5. CLEANUP & FORMATTING ---
    ds_out = xr.Dataset({
        'u_star': u_star,
        'tau': tau
    })

    ds_out.attrs['description'] = 'Bottom shear stress (tau) and friction velocity (u_star)'
    ds_out.attrs['model'] = 'SalishSeaCast'

    # --- 6. COMPUTE & SAVE (With Compression) ---
    output_file = '/ocean/vvalenzuela/MOAD/analysis-vicente/Ocean_Parcels/PBDE_simulations/bottom_shear_stress_2018.nc'
    
    # Tell NetCDF writer to compress the float32 arrays
    encoding_dict = {
        'u_star': {'zlib': True, 'complevel': 4, 'dtype': 'float32'},
        'tau':    {'zlib': True, 'complevel': 4, 'dtype': 'float32'}
    }

    print(f"Computing and streaming data to {output_file}...")

    # Optional: limit dask to smaller chunk sizes if it still peaks
    with dask.config.set({"array.slicing.split_large_chunks": True}):
        with ProgressBar():
            ds_out.to_netcdf(output_file, encoding=encoding_dict, engine='netcdf4')
        
    print("Computation complete. File saved!")

if __name__ == '__main__':
    main()