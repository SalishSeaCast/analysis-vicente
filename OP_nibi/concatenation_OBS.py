import xarray as xr
import numpy as np
from pathlib import Path

BASE_DIR = Path('/home/vicentev/scratch/vicentev/Simulation_V4/')

files = ['/home/vicentev/scratch/vicentev/Simulation_V4/Simulation_V4_year_1.nc',
         '/home/vicentev/scratch/vicentev/Simulation_V4/RESTARTS/R_2008/Simulation_V4_year_2.nc'
]

# Full paths
file_paths = [BASE_DIR / f for f in files]

# Concatenate along obs
ds = xr.open_mfdataset(
    file_paths,
    combine='nested',
    concat_dim='obs',
    )
ds = ds.assign_coords(obs=np.arange(len(ds.obs)))
# Save result
output_path = BASE_DIR / 'Simulation_V4_years_1_and_2.nc'
ds.to_netcdf(output_path)

print(f"Saved to {output_path}")