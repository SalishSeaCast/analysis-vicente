import xarray as xr
import numpy as np
from pathlib import Path

BASE_DIR = Path('/home/vicentev/scratch/vicentev/Simulation_V2/')

files = ['Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_V2_2008_P1.nc',
         'Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_V2_2008_P2.nc',
         'Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_V2_2008_P3.nc',
         'Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_V2_2008_P4.nc',
         'Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_V2_2008_P5.nc',
         'Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_V2_2008_P6.nc',
         'Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_V2_2008_P7.nc',
         'Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_V2_2008_P8.nc',
         'Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_V2_2008_P9.nc',
         'Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_V2_2008_P10.nc',
         'Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_V2_2008_P11.nc',
         'Tau_0_0025_Ads_0_05_MC_0_2_Vel_N_V2_2008_P12.nc',
]

# Full paths
file_paths = [BASE_DIR / f for f in files]

# Concatenate along obs
ds = xr.open_mfdataset(
    file_paths,
    combine='nested',
    concat_dim='trajectory',
    )
#
ds = ds.assign_coords(trajectory=np.arange(len(ds.trajectory)))
# Save result
output_path = BASE_DIR / 'Simulation_V2_year_2.nc'
ds.to_netcdf(output_path)

print(f"Saved to {output_path}")