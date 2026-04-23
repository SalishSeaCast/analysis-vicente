import xarray as xr
from pathlib import Path

BASE_DIR = Path('/home/vicentev/scratch/vicentev/Simulation_V3/')

files = [
    'PBDEs_0112007_run_365_days_Tau_0_005_Ads_0_05_MC_0_2_Vel_Hx1_2_P1.nc',
    'PBDEs_0112007_run_365_days_Tau_0_005_Ads_0_05_MC_0_2_Vel_Hx1_2_P2.nc',
    'PBDEs_0112007_run_365_days_Tau_0_005_Ads_0_05_MC_0_2_Vel_Hx1_2_P3.nc',
    'PBDEs_0112007_run_365_days_Tau_0_005_Ads_0_05_MC_0_2_Vel_Hx1_2_P4.nc',
    'PBDEs_0112007_run_365_days_Tau_0_005_Ads_0_05_MC_0_2_Vel_Hx1_2_P5.nc',
    'PBDEs_0112007_run_365_days_Tau_0_005_Ads_0_05_MC_0_2_Vel_Hx1_2_P6.nc',
    'PBDEs_0112007_run_365_days_Tau_0_005_Ads_0_05_MC_0_2_Vel_Hx1_2_P7.nc',
    'PBDEs_0112007_run_365_days_Tau_0_005_Ads_0_05_MC_0_2_Vel_Hx1_2_P8.nc',
    'PBDEs_0112007_run_365_days_Tau_0_005_Ads_0_05_MC_0_2_Vel_Hx1_2_P9.nc',
    'PBDEs_0112007_run_365_days_Tau_0_005_Ads_0_05_MC_0_2_Vel_Hx1_2_P10.nc',
    'PBDEs_0112007_run_365_days_Tau_0_005_Ads_0_05_MC_0_2_Vel_Hx1_2_P11.nc',
    'PBDEs_0112007_run_365_days_Tau_0_005_Ads_0_05_MC_0_2_Vel_Hx1_2_P12.nc',
]

# Full paths
file_paths = [BASE_DIR / f for f in files]

# Concatenate along obs
ds = xr.open_mfdataset(
    file_paths,
    combine='nested',
    concat_dim='obs'
)

# Save result
output_path = BASE_DIR / 'Simulation_V3_year_1.nc'
ds.to_netcdf(output_path)

print(f"Saved to {output_path}")