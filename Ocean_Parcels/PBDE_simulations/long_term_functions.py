# Functions to build long term times series and others
from matplotlib.path import Path as mPath
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon
from shapely.plotting import plot_polygon
import xarray as xr
from scipy.optimize import curve_fit
from matplotlib.path import Path
#
path = {'coords' : '/ocean/vvalenzuela/MOAD/grid/coordinates_seagrid_SalishSea201702.nc',
'mask' : '/ocean/vvalenzuela/MOAD/grid2/mesh_mask202108_TDV.nc',
'bathy' : '/ocean/vvalenzuela/MOAD/grid/bathymetry_202108.nc'
}
#
coords = xr.open_dataset(path['coords'], decode_times=False)
mask = xr.open_dataset(path['mask'])
bathy = xr.open_dataset(path['bathy'])
#
file_year1 = '/home/vvalenzuela/MOAD/Ocean_Parcels/Results_Final_Versions/Version_4_final/Simulation_V4_year_1.nc'
file_year2 = '/home/vvalenzuela/MOAD/Ocean_Parcels/Results_Final_Versions/Version_4_final/Simulation_V4_year_2.nc'
file_year3 = '/home/vvalenzuela/MOAD/Ocean_Parcels/Results_Final_Versions/Version_4_final/Simulation_V4_year_3.nc'
file_year4 = '/home/vvalenzuela/MOAD/Ocean_Parcels/Results_Final_Versions/Version_4_final/Simulation_V4_year_4.nc'
file_year5 = '/home/vvalenzuela/MOAD/Ocean_Parcels/Results_Final_Versions/Version_4_final/Simulation_V4_year_5.nc'
file_year6 = '/home/vvalenzuela/MOAD/Ocean_Parcels/Results_Final_Versions/Version_4_final/Simulation_V4_year_6.nc'
file_year7 = '/home/vvalenzuela/MOAD/Ocean_Parcels/Results_Final_Versions/Version_4_final/Simulation_V4_year_7.nc'
file_year8 = '/home/vvalenzuela/MOAD/Ocean_Parcels/Results_Final_Versions/Version_4_final/Simulation_V4_year_8.nc'
file_year9 = '/home/vvalenzuela/MOAD/Ocean_Parcels/Results_Final_Versions/Version_4_final/Simulation_V4_year_9.nc'
file_year10 = '/home/vvalenzuela/MOAD/Ocean_Parcels/Results_Final_Versions/Version_4_final/Simulation_V4_year_10.nc'
file_year11 = '/home/vvalenzuela/MOAD/Ocean_Parcels/Results_Final_Versions/Version_4_final/Simulation_V4_year_11.nc'
file_year12 = '/home/vvalenzuela/MOAD/Ocean_Parcels/Results_Final_Versions/Version_4_final/Simulation_V4_year_12.nc'
file_year13 = '/home/vvalenzuela/MOAD/Ocean_Parcels/Results_Final_Versions/Version_4_final/Simulation_V4_year_13.nc'
#
year1 = xr.open_dataset(file_year1)
year2 = xr.open_dataset(file_year2)
year3 = xr.open_dataset(file_year3)
year4 = xr.open_dataset(file_year4)
year5 = xr.open_dataset(file_year5)
year6 = xr.open_dataset(file_year6)
year7 = xr.open_dataset(file_year7)
year8 = xr.open_dataset(file_year8)
year9 = xr.open_dataset(file_year9)
year10 = xr.open_dataset(file_year10)
year11 = xr.open_dataset(file_year11)
year12 = xr.open_dataset(file_year12)
year13 = xr.open_dataset(file_year13)
#
def pad_year(year, conversion, ts_full_yr1):
    ts_full_yeard = (year.status == 2).sum(axis=0) * conversion
    ts_full_year = np.zeros_like(ts_full_yr1)
    ts_full_year[1::2] = ts_full_yeard
    ts_full_year[0] = 0.5 * (ts_full_year[1] + ts_full_year[-1])
    ts_full_year[2::2] = 0.5 * (ts_full_year[1:-2:2] + ts_full_year[3::2])
    return ts_full_year
#
conversion = 1 / 96 # 1 g per day = 96 particles
ts_full_yr1 = (year1.status == 2).sum(axis=0) * conversion
ts_full_yr2 = (year2.status == 2).sum(axis=0) * conversion
ts_full_yr3 = pad_year(year3, conversion, ts_full_yr1)
ts_full_yr4 = pad_year(year4, conversion, ts_full_yr1)
ts_full_yr5 = pad_year(year5, conversion, ts_full_yr1)
ts_full_yr6 = pad_year(year6, conversion, ts_full_yr1)
ts_full_yr7 = pad_year(year7, conversion, ts_full_yr1)
ts_full_yr8 = pad_year(year8, conversion, ts_full_yr1)
ts_full_yr9 = pad_year(year9, conversion, ts_full_yr1)
ts_full_yr10 = pad_year(year10, conversion, ts_full_yr1)
ts_full_yr11 = pad_year(year11, conversion, ts_full_yr1)
ts_full_yr12 = pad_year(year12, conversion, ts_full_yr1)
ts_full_yr13 = pad_year(year13, conversion, ts_full_yr1)
#
max_year = 2050
#
def total_pool_tseries(congener_min, congener_max, BDE_name):
    #
    BDE = pd.DataFrame({
    'Date': pd.to_datetime(congener_min['Unnamed: 0']),
    'Min': congener_min[BDE_name],
    'Max': congener_max[BDE_name]
    })
    #
    BDE['year'] = BDE['Date'].dt.year
    discharge_min = BDE.groupby('year')['Min'].mean()
    discharge_max = BDE.groupby('year')['Max'].mean()
    # MAXIMUM
    one_year = 365 * 8
    zero_values = np.zeros((max_year-1970+1)*one_year)
    one_values = np.zeros_like(zero_values)
    two_values = np.zeros_like(zero_values)
    three_values = np.zeros_like(zero_values)
    four_values = np.zeros_like(zero_values)
    five_values = np.zeros_like(zero_values)
    six_values = np.zeros_like(zero_values)
    seven_values = np.zeros_like(zero_values)
    eight_values = np.zeros_like(zero_values)
    nine_values = np.zeros_like(zero_values)
    ten_values = np.zeros_like(zero_values)
    eleven_values = np.zeros_like(zero_values)
    twelve_values = np.zeros_like(zero_values)
    #
    for year in range(1970, max_year):
        zero_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            discharge_max[year] * ts_full_yr1)
    for year in range(1971, max_year):
        one_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            zero_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-1] * ts_full_yr2)
    for year in range(1972, max_year):
        two_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            one_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-2] * ts_full_yr3)
    for year in range(1973, max_year):
        three_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            two_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-3] * ts_full_yr4)
    for year in range(1974, max_year):
        four_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            three_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-4] * ts_full_yr5)
    for year in range(1975, max_year):
        five_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            four_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-5] * ts_full_yr6)
    for year in range(1976, max_year):
        six_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            five_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-6] * ts_full_yr7)
    for year in range(1977, max_year):
        seven_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            six_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-7] * ts_full_yr8)
    for year in range(1978, max_year):
        eight_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            seven_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-8] * ts_full_yr9)    
    for year in range(1979, max_year):
        nine_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            eight_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-9] * ts_full_yr10)
    for year in range(1980, max_year):
        ten_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            nine_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-10] * ts_full_yr11)           
    for year in range(1981, max_year):
        eleven_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            ten_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-11] * ts_full_yr12)
    for year in range(1982, max_year):
        twelve_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            eleven_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-12] * ts_full_yr13)
    #
    years_tseries_max = [zero_values, one_values, two_values, three_values, four_values, five_values,
                         six_values, seven_values, eight_values, nine_values, ten_values, eleven_values, twelve_values]
    # MINIMUM     
    zero_values = np.zeros((max_year-1970+1)*one_year)
    one_values = np.zeros_like(zero_values)
    two_values = np.zeros_like(zero_values)
    three_values = np.zeros_like(zero_values)
    four_values = np.zeros_like(zero_values)
    five_values = np.zeros_like(zero_values)
    six_values = np.zeros_like(zero_values)
    seven_values = np.zeros_like(zero_values)
    eight_values = np.zeros_like(zero_values)
    nine_values = np.zeros_like(zero_values)
    ten_values = np.zeros_like(zero_values)
    eleven_values = np.zeros_like(zero_values)
    twelve_values = np.zeros_like(zero_values)
    #
    for year in range(1970, max_year):
        zero_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            discharge_min[year] * ts_full_yr1)
    for year in range(1971, max_year):
        one_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            zero_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-1] * ts_full_yr2)
    for year in range(1972, max_year):
        two_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            one_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-2] * ts_full_yr3)
    for year in range(1973, max_year):
        three_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            two_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-3] * ts_full_yr4)
    for year in range(1974, max_year):
        four_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            three_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-4] * ts_full_yr5)
    for year in range(1975, max_year):
        five_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            four_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-5] * ts_full_yr6)
    for year in range(1976, max_year):
        six_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            five_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-6] * ts_full_yr7)
    for year in range(1977, max_year):
        seven_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            six_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-7] * ts_full_yr8)
    for year in range(1978, max_year):
        eight_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            seven_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-8] * ts_full_yr9)    
    for year in range(1979, max_year):
        nine_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            eight_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-9] * ts_full_yr10)
    for year in range(1980, max_year):
        ten_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            nine_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-10] * ts_full_yr11)           
    for year in range(1981, max_year):
        eleven_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            ten_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-11] * ts_full_yr12)
    for year in range(1982, max_year):
        twelve_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            eleven_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-12] * ts_full_yr13)
    #
    years_tseries_min = [zero_values, one_values, two_values, three_values, four_values, five_values,
                         six_values, seven_values, eight_values, nine_values, ten_values, eleven_values, twelve_values]                
    #
    return years_tseries_min, years_tseries_max             
#
# Time series per region
def get_total_timeseries_optimized(filename, polygon_coords, 
                                  lon_var='lon', lat_var='lat', status_var='status', 
                                  target_status=2, chunk_size=500):
    """
    Computes a particle timeseries filtered by status and polygon containment.
    Optimized for minimal memory usage and high execution speed via chunking.
    """
    poly_path = Path(polygon_coords)
    
    with xr.open_dataset(filename) as data:
        n_obs = data.sizes['obs'] 
        counts_per_time = np.zeros(n_obs, dtype=int)
        
        for start in range(0, n_obs, chunk_size):
            end = min(start + chunk_size, n_obs)
            
            status_chunk = data[status_var].isel(obs=slice(start, end)).values
            
            if isinstance(target_status, (list, tuple, np.ndarray)):
                status_mask = np.isin(status_chunk, target_status)
            else:
                status_mask = (status_chunk == target_status)
            
            if not np.any(status_mask):
                continue
                

            lon_chunk = data[lon_var].isel(obs=slice(start, end)).values[status_mask]
            lat_chunk = data[lat_var].isel(obs=slice(start, end)).values[status_mask]
            
            points = np.column_stack((lon_chunk, lat_chunk))
            inside_poly = poly_path.contains_points(points)
            
            if np.any(inside_poly):

                chunk_time_indices = np.where(status_mask)[1]
                absolute_time_indices = chunk_time_indices[inside_poly] + start
                
                counts_per_time += np.bincount(absolute_time_indices, minlength=n_obs)
                
    return pd.DataFrame({'Total_Count': counts_per_time})
#
#
def regions_tseries(congener_min, congener_max, polygon_lon_lat, BDE_name):
    #
    t_series_region = [get_total_timeseries_optimized(file_year1, polygon_lon_lat)['Total_Count'].values[::2],
                       get_total_timeseries_optimized(file_year2, polygon_lon_lat)['Total_Count'].values[::2],
                       get_total_timeseries_optimized(file_year3, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year4, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year5, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year6, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year7, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year8, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year9, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year10, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year11, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year12, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year13, polygon_lon_lat)['Total_Count'].values
                       ]
    #
    BDE = pd.DataFrame({
    'Date': pd.to_datetime(congener_min['Unnamed: 0']),
    'Min': congener_min[BDE_name],
    'Max': congener_max[BDE_name]
    })
    #
    BDE['year'] = BDE['Date'].dt.year
    discharge_min = BDE.groupby('year')['Min'].mean()
    discharge_max = BDE.groupby('year')['Max'].mean()
    # MAXIMUM
    one_year = 365 * 4
    zero_values = np.zeros((max_year-1970+1)*one_year)
    one_values = np.zeros_like(zero_values)
    two_values = np.zeros_like(zero_values)
    three_values = np.zeros_like(zero_values)
    four_values = np.zeros_like(zero_values)
    five_values = np.zeros_like(zero_values)
    six_values = np.zeros_like(zero_values)
    seven_values = np.zeros_like(zero_values)
    eight_values = np.zeros_like(zero_values)
    nine_values = np.zeros_like(zero_values)
    ten_values = np.zeros_like(zero_values)
    eleven_values = np.zeros_like(zero_values)
    twelve_values = np.zeros_like(zero_values)
    #
    for year in range(1970, max_year):
        zero_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            discharge_max[year] * t_series_region[0])
    for year in range(1971, max_year):
        one_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            zero_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-1] * t_series_region[1])
    for year in range(1972, max_year):
        two_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            one_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-2] * t_series_region[2])
    for year in range(1973, max_year):
        three_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            two_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-3] * t_series_region[3])
    for year in range(1974, max_year):
        four_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            three_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-4] * t_series_region[4])
    for year in range(1975, max_year):
        five_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            four_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-5] * t_series_region[5])
    for year in range(1976, max_year):
        six_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            five_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-6] * t_series_region[6])
    for year in range(1977, max_year):
        seven_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            six_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-7] * t_series_region[7])
    for year in range(1978, max_year):
        eight_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            seven_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-8] * t_series_region[8])    
    for year in range(1979, max_year):
        nine_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            eight_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-9] * t_series_region[9])
    for year in range(1980, max_year):
        ten_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            nine_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-10] * t_series_region[10])           
    for year in range(1981, max_year):
        eleven_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            ten_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-11] * t_series_region[11])
    for year in range(1982, max_year):
        twelve_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            eleven_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_max[year-12] * t_series_region[12])
    #
    years_tseries_max = [zero_values, one_values, two_values, three_values, four_values, five_values,
                         six_values, seven_values, eight_values, nine_values, ten_values, eleven_values, twelve_values]
    # MINIMUM     
    zero_values = np.zeros((max_year-1970+1)*one_year)
    one_values = np.zeros_like(zero_values)
    two_values = np.zeros_like(zero_values)
    three_values = np.zeros_like(zero_values)
    four_values = np.zeros_like(zero_values)
    five_values = np.zeros_like(zero_values)
    six_values = np.zeros_like(zero_values)
    seven_values = np.zeros_like(zero_values)
    eight_values = np.zeros_like(zero_values)
    nine_values = np.zeros_like(zero_values)
    ten_values = np.zeros_like(zero_values)
    eleven_values = np.zeros_like(zero_values)
    twelve_values = np.zeros_like(zero_values)
    #
    for year in range(1970, max_year):
        zero_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            discharge_min[year] * t_series_region[0])
    for year in range(1971, max_year):
        one_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            zero_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-1] * t_series_region[1])
    for year in range(1972, max_year):
        two_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            one_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-2] * t_series_region[2])
    for year in range(1973, max_year):
        three_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            two_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-3] * t_series_region[3])
    for year in range(1974, max_year):
        four_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            three_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-4] * t_series_region[4])
    for year in range(1975, max_year):
        five_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            four_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-5] * t_series_region[5])
    for year in range(1976, max_year):
        six_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            five_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-6] * t_series_region[6])
    for year in range(1977, max_year):
        seven_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            six_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-7] * t_series_region[7])
    for year in range(1978, max_year):
        eight_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            seven_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-8] * t_series_region[8])    
    for year in range(1979, max_year):
        nine_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            eight_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-9] * t_series_region[9])
    for year in range(1980, max_year):
        ten_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            nine_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-10] * t_series_region[10])           
    for year in range(1981, max_year):
        eleven_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            ten_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-11] * t_series_region[11])
    for year in range(1982, max_year):
        twelve_values[(year-1970)*one_year:(year-1970+1)*one_year] = (
            eleven_values[(year-1970)*one_year:(year-1970+1)*one_year] + 
            discharge_min[year-12] * t_series_region[12])
    #
    years_tseries_min = [zero_values, one_values, two_values, three_values, four_values, five_values,
                         six_values, seven_values, eight_values, nine_values, ten_values, eleven_values, twelve_values]                
    #
    return years_tseries_min, years_tseries_max
#
# EXTRAPOLATION FUNCTIONS
#
def curve_double_yr(t, maxv1, maxv2, decay1, decay2, tmm=1):
    value = maxv1 * np.exp(-t/decay1) * (np.exp(np.minimum(t, tmm)/decay1) - 1)
    value = value + maxv2 * np.exp(-t/decay2) * (np.exp(np.minimum(t, tmm)/decay2) - 1)
    return value
#
# Extrapolation Total SoG
def extrapolation_total(congener_min, congener_max, BDE_name):
    tseries_concat = np.concatenate([ts_full_yr1, ts_full_yr2, ts_full_yr3, ts_full_yr4, ts_full_yr5,
                                    ts_full_yr6, ts_full_yr7, ts_full_yr8, ts_full_yr9, ts_full_yr10,
                                    ts_full_yr11, ts_full_yr12, ts_full_yr13]) * conversion
    #
    fit_double_yr, cov = curve_fit(curve_double_yr, np.arange(0, 13, 1/2920), tseries_concat, p0=[0.6, 0.3, 2, 5])
    #
    BDE = pd.DataFrame({
    'Date': pd.to_datetime(congener_min['Unnamed: 0']),
    'Min': congener_min[BDE_name],
    'Max': congener_max[BDE_name]
    })
    #
    BDE['year'] = BDE['Date'].dt.year
    discharge_min = BDE.groupby('year')['Min'].mean()
    discharge_max = BDE.groupby('year')['Max'].mean()
    #
    # Min
    pred_min = discharge_min[1970] * curve_double_yr(np.arange(0, max_year-1970, 1/2920), 
                                          fit_double_yr[0], fit_double_yr[1], fit_double_yr[2], fit_double_yr[3],)
    for year in range(1970+1, max_year):
        pred_min[2920*(year-1970):] = (pred_min[2920*(year-1970):] + 
                discharge_min[year] * curve_double_yr(np.arange(0, max_year-year, 1/2920), 
                                                fit_double_yr[0], fit_double_yr[1], fit_double_yr[2], fit_double_yr[3]))
    # Max
    pred_max = discharge_max[1970] * curve_double_yr(np.arange(0, max_year-1970, 1/2920), 
                                          fit_double_yr[0], fit_double_yr[1], fit_double_yr[2], fit_double_yr[3],)
    for year in range(1970+1, max_year):
        pred_max[2920*(year-1970):] = (pred_max[2920*(year-1970):] + 
                discharge_max[year] * curve_double_yr(np.arange(0, max_year-year, 1/2920), 
                                                fit_double_yr[0], fit_double_yr[1], fit_double_yr[2], fit_double_yr[3])) 
    return pred_min, pred_max        
#
# Regional Extrapolation
#     
def extrapolation_regional(congener_min, congener_max, polygon_lon_lat, BDE_name):
    #
    tseries = [get_total_timeseries_optimized(file_year1, polygon_lon_lat)['Total_Count'].values[::2],
                       get_total_timeseries_optimized(file_year2, polygon_lon_lat)['Total_Count'].values[::2],
                       get_total_timeseries_optimized(file_year3, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year4, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year5, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year6, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year7, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year8, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year9, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year10, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year11, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year12, polygon_lon_lat)['Total_Count'].values,
                       get_total_timeseries_optimized(file_year13, polygon_lon_lat)['Total_Count'].values
                       ]
    #
    tseries_region = np.concatenate([tseries[0], tseries[1], tseries[2], tseries[3], tseries[4], tseries[5], tseries[6]
                                     , tseries[7], tseries[8], tseries[9], tseries[10], tseries[11], tseries[12]])   
    #
    fit_double_yr, cov = curve_fit(curve_double_yr, np.arange(0, 13, 1/1460), tseries_region, p0=[0.6, 0.3, 2, 5])
    #
    BDE = pd.DataFrame({
    'Date': pd.to_datetime(congener_min['Unnamed: 0']),
    'Min': congener_min[BDE_name],
    'Max': congener_max[BDE_name]
    })
    #
    BDE['year'] = BDE['Date'].dt.year
    discharge_min = BDE.groupby('year')['Min'].mean()
    discharge_max = BDE.groupby('year')['Max'].mean()

    #
    # Min
    pred_min = discharge_min[1970] * curve_double_yr(np.arange(0, max_year-1970, 1/1460), 
                                          fit_double_yr[0], fit_double_yr[1], fit_double_yr[2], fit_double_yr[3],)
    for year in range(1970+1, max_year):
        pred_min[1460*(year-1970):] = (pred_min[1460*(year-1970):] + 
                discharge_min[year] * curve_double_yr(np.arange(0, max_year-year, 1/1460), 
                                                fit_double_yr[0], fit_double_yr[1], fit_double_yr[2], fit_double_yr[3]))
    # Max
    pred_max = discharge_max[1970] * curve_double_yr(np.arange(0, max_year-1970, 1/1460), 
                                          fit_double_yr[0], fit_double_yr[1], fit_double_yr[2], fit_double_yr[3],)
    for year in range(1970+1, max_year):
        pred_max[1460*(year-1970):] = (pred_max[1460*(year-1970):] + 
                discharge_max[year] * curve_double_yr(np.arange(0, max_year-year, 1/1460), 
                                                fit_double_yr[0], fit_double_yr[1], fit_double_yr[2], fit_double_yr[3])) 
    return pred_min, pred_max      
