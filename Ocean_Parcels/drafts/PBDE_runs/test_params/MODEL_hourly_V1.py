#
# PBDEs Dispersion Model based on Ocean Parcels and Salish Sea Cast
#
# This model describes PBDEs dynamics including physical and chemical properties.
#
# An initial source point is located in the Iona Outfall diffusers coordinates at around 70 m depth.
#
# This version releases particles every hour for a certain amount of days (INPUT)
#
import sys
import random
import xarray as xr
import numpy as np
import os
import yaml
import math
from datetime import datetime, timedelta
from parcels import FieldSet, Field, VectorField, ParticleSet, JITParticle, ParcelsRandom, Variable, Kernel, AdvectionRK4
from cartopy import crs, feature
import zarr 
#
sys.path.append('/ocean/vvalenzuela/MOAD/Ocean_Parcels')
#
from OP_functions import *
#
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
#
#
def PBDEs_OP_run(start_time, sim_length, number_particles, days_release, delta_t):
    #### DATA AND OUTPUT PATHS ####
    path = {'NEMO': '/results2/SalishSea/nowcast-green.202111/',
    'coords': '/ocean/vvalenzuela/MOAD/grid/coordinates_seagrid_SalishSea201702.nc',
    'coordsWW3': '/ocean/vvalenzuela/MOAD/grid2/WW3_grid.nc',
    'mask': '/ocean/vvalenzuela/MOAD/grid2/mesh_mask202108_TDV.nc',
    'bat': '/ocean/vvalenzuela/MOAD/grid/bathymetry_202108.nc',
    'out': '/home/vvalenzuela/MOAD/Ocean_Parcels/results/Test_runs/Test_hourly_V2',
    'home': '/home/vvalenzuela/MOAD/Ocean_Parcels',
    'anim': '/home/vvalenzuela/MOAD/Ocean_Parcels/results/PBDE_runs/animations'}
    #
    coords = xr.open_dataset(path['coords'], decode_times=False)
    mask = xr.open_dataset(path['mask'])
    #### Function for getting dates stamps ####
    def get_timestamps(start_time,sim_length):
        timestamps=[]
        duration = timedelta(days=sim_length)
        for day in range(duration.days):
            timestamps.append([start_time + timedelta(days=day)])
        return np.array(timestamps, dtype='datetime64')
    #
    #### Function to get grid point gridX and gridY ####
    path_NEMO = make_prefix(start_time,path['NEMO'])
    jjii = xr.open_dataset('/ocean/vvalenzuela/MOAD/grid/grid_from_lat_lon_mask999.nc')
    def finder(lati,loni):
        j = [jjii.jj.sel(lats=lati, lons=loni, method='nearest').item()][0]
        i = [jjii.ii.sel(lats=lati, lons=loni, method='nearest').item()][0]
        return j,i
    #
    #### Setting deploying cordinates ####
    #
    # IONA OUTFALL COORDINATES
    clat = [49.195045]
    clon = [-123.301956]
    z = 70 # 70 m depth
    #
    n_hourly = int(number_particles) # amount of particles released every hour through the day
    #
    release_interval = timedelta(hours=1)
    total_release_duration = timedelta(days=days_release)
    #
    release_times = [start_time + i * release_interval for i in range(int(total_release_duration / release_interval))]
    num_releases = len(release_times)
    #
    lon = [clon[0]] * num_releases * n_hourly
    lat = [clat[0]] * num_releases * n_hourly
    depth = [z] * num_releases * n_hourly
    time = [release_time for release_time in release_times for _ in range(n_hourly)]
    #
    a, b = finder(clat[0], clon[0])
    #print ("The total depth at this location is", mask.totaldepth[a, b].values, 'm')
    #
    duration = timedelta(days=sim_length) # RUN DURATION IN DAYS
    #
    #
    #### Name of the output file #### 
    name_states = 'PBDEs_run_for_'+str(sim_length)+'_days_'+str(n_hourly)+'_hourly_particles_'
    daterange = [start_time+timedelta(days=i) for i in range(sim_length)]
    fn =  name_states + '_'.join(d.strftime('%Y%m%d')+'_1n' for d in [start_time, start_time+duration]) + '.zarr'
    outfile_states = os.path.join(path['out'], fn)
    #
    local = 0
    ####
    ####
    ####
    #### CREATING FIELDSETS ####
    varlist=['U','V','W']
    filenames,variables=filename_set(start_time,sim_length,varlist)
    dimensions = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw','time': 'time_counter'}
    field_set=FieldSet.from_nemo(filenames, variables, dimensions, allow_time_extrapolation=True, chunksize='auto')
    #
    #Find file names and variable names ###'Diat','Flag'###
    varlist=['US','VS','WL','R','T','S','ssh','Bathy','Kz','totdepth','Vol','last_cell_index']
    filenames,variables=filename_set(start_time,sim_length,varlist)
    dimensions = {'lon': 'glamt', 'lat': 'gphit', 'depth': 'deptht', 'time': 'time_counter'}
    density = Field.from_netcdf(filenames['R'], variables['R'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    field_set.add_field(density)
    #
    #Add Vertical diffusivity coefficient field
    dimensions = {'lon': 'glamt', 'lat': 'gphit', 'depth': 'depthw','time': 'time_counter'}
    Kz = Field.from_netcdf(filenames['Kz'], variables['Kz'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    field_set.add_field(Kz)
    #
    #Add Bathymetry 2D field
    dimensions = {'lon': 'glamt', 'lat': 'gphit'}
    Bth = Field.from_netcdf(filenames['Bathy'], variables['Bathy'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    TD = Field.from_netcdf(filenames['totdepth'], variables['totdepth'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    MBATHY = Field.from_netcdf(filenames['last_cell_index'], variables['last_cell_index'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    field_set.add_field(MBATHY)
    field_set.add_field(Bth)
    field_set.add_field(TD)
    #
    #Add SSH 
    dimensions = {'lon': 'glamt', 'lat': 'gphit','time': 'time_counter'}
    SSH = Field.from_netcdf(filenames['ssh'], variables['ssh'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    field_set.add_field(SSH)
    #
    # Add e3t
    varlist = ['cell_size']
    filenames,variables=filename_set(start_time,sim_length,varlist)
    dimensions = {'lon': 'glamt', 'lat': 'gphit', 'depth': 'deptht', 'time': 'time_counter'}
    E3T = Field.from_netcdf(filenames['cell_size'], variables['cell_size'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    field_set.add_field(E3T)
    ####
    #### DEFINE A PARTICLE TYPE AND SET ####
    class MPParticle(JITParticle):    
        n = Variable('n', initial = n_hourly)
        vvl_factor = Variable('fact', initial =  1)    
        wa = Variable('wa', initial =  0) 
        wm = Variable('wm', initial =  0)
        initialized = Variable('initialized', initial = 0)
        status = Variable('status') # different status for different processes
    #
    pset_states = ParticleSet.from_list(field_set, MPParticle, lon=lon, lat=lat,
    depth=depth, time=time)           
    ##########################################################
    #
    #
    runtime = duration  # Total simulation runtime
    dt = delta_t     # Simulation timestep in seconds
    output_interval = timedelta(hours=1)
    output_file = pset_states.ParticleFile(name=outfile_states, outputdt=output_interval)
    kernels = pset_states.Kernel(PBDEs_states) +  pset_states.Kernel(PBDEs_forms) + pset_states.Kernel(Advection) + pset_states.Kernel(turb_mix) + pset_states.Kernel(Displacement) + pset_states.Kernel(resuspension) + pset_states.Kernel(CheckOutOfBounds) + pset_states.Kernel(KeepInOcean)
    #
    #
    pset_states.execute(kernels,
            runtime=runtime,
            dt=dt,
            output_file=output_file)
    #
    print('Output file at ', outfile_states)
        
if __name__ == "__main__":
    print("RUNNING PBDEs_OP_run! :D")
    # Input from the terminal
    start_time_str = sys.argv[1]  # Example: "2022-01-01T00:00:00"
    length_sim = int(sys.argv[2]) # Example: 10
    number_of_particles_per_release = int(sys.argv[3]) # Example: 5  
    number_of_days_release = int(sys.argv[4]) # Example: 2   
    delta_t = int(sys.argv[5]) # Example: 40 [in seconds]
    # Convert start_time_str to datetime
    start_time = datetime.strptime(start_time_str, "%Y-%m-%dT%H:%M:%S")
    #
    #### this is for avoiding time overlapping #### 
    total_days = length_sim + (number_of_days_release*2)
    PBDEs_OP_run(start_time, total_days, number_of_particles_per_release, number_of_days_release, delta_t)
    #
    print(f"PBDEs_OP_run simulation for {total_days} days, releasing {number_of_particles_per_release} particles every 1 hour for {number_of_days_release} days from the starting date {start_time}")
    #
    ## How to run in the terminal:
    # python -m MODEL_hourly_V1 start_time_str length_sim number_particles days_release delta_t