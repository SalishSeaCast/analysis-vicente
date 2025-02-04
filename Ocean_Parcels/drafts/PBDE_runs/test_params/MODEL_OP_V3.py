# RECONSTRUCTED PBDEs Dispersion Model VERSION v3
#
#
# THIS VERSION WORKS MUCH BETTER #
# Initialization set with particle.initialized = 0 ---> Works good
# dt_h = 1 / 3600 to be consistent with the seconds conversion for any particle.dt
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
####
####
####
####
start = datetime(2022, 1, 1) #Start date
length = 6 # Set Time length [days] 
dt = 1#time step in seconds 
N = 1 # number of  locations
n = 100 # number of particles per location
dmin = 60 #minimum depth
dd = 20 #max depth difference from dmin
dtp = 0
odt = 1 
rrr = 1e3
####
####
####
####
def simulation_run(start_time, sim_length, number_particles, dt_resolution):
    #
    #
    ####
    #### DATA AND OUTPUT PATHS ####
    path = {'NEMO': '/results2/SalishSea/nowcast-green.202111/',
    'coords': '/ocean/vvalenzuela/MOAD/grid/coordinates_seagrid_SalishSea201702.nc',
    'coordsWW3': '/ocean/vvalenzuela/MOAD/grid2/WW3_grid.nc',
    'mask': '/ocean/vvalenzuela/MOAD/grid2/mesh_mask202108_TDV.nc',
    'bat': '/ocean/vvalenzuela/MOAD/grid/bathymetry_202108.nc',
    'out': '/home/vvalenzuela/MOAD/Ocean_Parcels/results/Test_runs',
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
    #
    a, b = finder(clat[0], clon[0])
    print ("The total depth at this location is", mask.totaldepth[a, b].values, 'm')
    #
    duration = timedelta(days=sim_length) # RUN DURATION IN DAYS
    #
    x_offset, y_offset, z = p_deploy(N,number_particles,dmin,dd,rrr)
    #
    lon = np.zeros([N,number_particles])
    lat = np.zeros([N,number_particles])
    for i in range(N):
        lon[i,:]=(clon[i] + x_offset[i,:])
        lat[i,:]=(clat[i] + y_offset[i,:])
    #
    #
    #### Name of the output file #### 
    name_states = 'PBDEs_run_for_'+str(sim_length)+'_days_'+str(number_particles)+'_particles_'+str(dt_resolution)+'s_dt_' 
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
    ####
    ####
    #### DEFINE A PARTICLE TYPE AND SET ####
    class MPParticle(JITParticle):    
        n = Variable('n', initial = number_particles)
        vvl_factor = Variable('fact', initial =  1)    
        wa = Variable('wa', initial =  0) 
        wm = Variable('wm', initial =  0)
        initialized = Variable('initialized', initial = 0)
        status = Variable('status') # different status for different processes
    #
    pset_states = ParticleSet.from_list(field_set, MPParticle, lon=lon, lat=lat, depth=z, time=start_time+timedelta(hours=odt))
    ####
    ####
    ####
    ####
    ################ KERNELS ###################
    #
    #### States Kernel for particles release ####
    #### Interaction of 3 different states for PBDEs ####
    # 
    def PBDEs_states(particle, fieldset, time):
        if particle.initialized == 0:#particle.time <= 3600
            n = particle.n 
            # n is the total amount of particles released at the starting location
            data = ParcelsRandom.randint(0, n-1)
            #
            # PBDEs as Sewage Particles
            if data < 3*(n/4):
                particle.status = 1
            #
            # Colloidal/Dissolved PBDEs
            else:
                particle.status = 2
            #
            particle.initialized = 1
            #print('Particle Initialized = 0')    
        else: #particle.time > 3600:
            abso = 0.7/(24)#0.038/(24) #per hour
            deso_s = 3.2/(24) #per hour
            deso_m = 1.8/(24)#1.6/(24) #per hour
            #dt_h = 1 / 3600
            dt_h = (math.fabs(particle.dt)/math.fabs(particle.dt)) / 3600 # forcing to have a 1 second resolution
            value = ParcelsRandom.random()# * dt_h
            #
            #value = ParcelsRandom.random() * dt_h
            if particle.status == 2:
                if value < 1 - math.exp(-abso * dt_h):
                    particle.status = 3
               # From Coloidal/Dissolved form to being attached to a Marine Particle           
            elif particle.status == 1:
                if value < 1 - math.exp(-deso_s * dt_h):
                    particle.status = 2
               # From Sewage Particle to Colloidal/Dissolved PBDE form
            elif particle.status == 3:
                if value < 1 - math.exp(-deso_m * dt_h):
                    particle.status = 2
    #
    #### DYNAMICS KERNELS KEEP particle.dt as it is: If I forced to be always =1, then there will be inconsistencies between processes!!!
    #### PBDEs states sinking velocities features ####
    def PBDEs_forms(particle, fieldset, time):
        #
        #dt_h = 1 / 3600
        dt_h = (math.fabs(particle.dt)/math.fabs(particle.dt)) / 3600 # forcing to have a 1 second resolution
        if particle.status == 1:
            sinkvel = 50*(dt_h) # m/hr * dt --> to seconds
            particle.depth += sinkvel * particle.dt
        #Sewage Particles sink fast        
        elif particle.status == 2:
            sinkvel = 0.0
            particle.depth += sinkvel * particle.dt
        # Colloids just float around and move with advection
        elif particle.status == 3:
            sinkvel = 10*(dt_h) # m/hr * dt --> to seconds
            particle.depth += sinkvel * particle.dt 
    #
    #### ADVECTION ####
    def Advection(particle, fieldset, time): 
        # Advection for all PBDEs in status 1, 2 and 3
        if particle.status == 1 or particle.status == 2 or particle.status == 3: 
            ssh = fieldset.sossheig[time, particle.depth, particle.lat, particle.lon] #SSH(t) sea surface height
            sshn = fieldset.sossheig[time+particle.dt, particle.depth, particle.lat, particle.lon] #SSH(t+dt) sea surface height in the next time step
            td = fieldset.totaldepth[time, particle.depth, particle.lat, particle.lon]#Total_depth 
            particle.fact = (1+ssh/td)
            VVL = (sshn-ssh)*particle.depth/(td)
            #VVL = (sshn-ssh)*particle.depth/(td+ssh)
            (u1, v1, w1) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
            lon1 = particle.lon + u1*.5*particle.dt
            lat1 = particle.lat + v1*.5*particle.dt
            dep1 = particle.depth + w1*.5*particle.dt/particle.fact
            (u2, v2, w2) = fieldset.UVW[time + .5 * particle.dt, dep1, lat1, lon1]
            lon2 = particle.lon + u2*.5*particle.dt
            lat2 = particle.lat + v2*.5*particle.dt
            dep2 = particle.depth + w2*.5*particle.dt/particle.fact
            (u3, v3, w3) = fieldset.UVW[time + .5 * particle.dt, dep2, lat2, lon2]
            lon3 = particle.lon + u3*particle.dt
            lat3 = particle.lat + v3*particle.dt
            dep3 = particle.depth + w3*particle.dt/particle.fact
            (u4, v4, w4) = fieldset.UVW[time + particle.dt, dep3, lat3, lon3]
            wa = (w1 + 2*w2 + 2*w3 + w4) /6.
            particle.wa = wa* particle.dt
            particle_dlon = (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
            particle_dlat = (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
            particle_ddepth = particle.wa/particle.fact + VVL
            if particle_ddepth + particle.depth < 0:
                particle_ddepth = - (2*particle.depth+particle_ddepth)
#        else:
#            particle_dlon = 0
#            particle_dlat = 0
#            particle_ddepth = 0     
    #
    #### TURBULENT MIX ####
    def turb_mix(particle,fieldset,time):
        if particle.status == 1 or particle.status == 2 or particle.status == 3:
            """Vertical mixing"""
            #Vertical mixing
            if particle.depth + 0.5/particle.fact > td: #Only calculate gradient of diffusion for particles deeper than 0.5 otherwise OP will check for particles outside the domain and remove it.
                Kzdz = 2*(fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon]-fieldset.vert_eddy_diff[time, particle.depth-0.5/particle.fact, particle.lat, particle.lon]) #backwards difference 
            else: 
                Kzdz = 2*(fieldset.vert_eddy_diff[time, particle.depth+0.5/particle.fact, particle.lat, particle.lon]-fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon]) #forward difference 
            dgrad = Kzdz*particle.dt/particle.fact
            if particle.depth+(0.5*dgrad) > 0 and particle.depth+(0.5*dgrad) < td:
                Kz = fieldset.vert_eddy_diff[time, particle.depth+ 0.5*dgrad, particle.lat, particle.lon] #Vertical diffusivity SSC  
            else:
                Kz = 0 
            #
            Rr = ParcelsRandom.uniform(-1, 1)
            d_random = sqrt(3*2*Kz*particle.dt) * Rr/particle.fact
            dzs = (dgrad + d_random)
            particle.wm = dzs*particle.fact
    #
    #### VERTICAL DISPLACEMENT ####
    def Displacement(particle,fieldset,time):
        if particle.status == 1 or particle.status == 2 or particle.status == 3:
            #Apply turbulent mixing.
            if dzs + particle_ddepth + particle.depth > td:
                particle.depth  = td # Get particles attached to the bottom when they reach it
                particle.status = 4
            #
            elif dzs + particle.depth+ particle_ddepth < 0:
                particle_ddepth = -(dzs + particle.depth+particle_ddepth) #reflection on surface
            #
            else:
                particle_ddepth += dzs #apply mixing
    #
    #### RESUSPENSION ####
    def resuspension(particle, fieldset, time):
        if particle.status == 4:
            threshold = 1 # threshold for particles to know when to resuspend
            # Calculation of U_star, which is proportional to the bottom stress (tau)
            k = 0.42
            z_star = 0.07
            u_horizontal = (1/4) * (fieldset.U[time, fieldset.mbathy[time, particle.depth, particle.lat, particle.lon] - 1, particle.lat, particle.lon] + fieldset.U[time, fieldset.mbathy[time, particle.depth, particle.lat, particle.lon] - 1, particle.lat, particle.lon -1]) ** 2
            v_horizontal = (1/4) * (fieldset.V[time, fieldset.mbathy[time, particle.depth, particle.lat, particle.lon] - 1, particle.lat, particle.lon] + fieldset.U[time, fieldset.mbathy[time, particle.depth, particle.lat, particle.lon] - 1, particle.lat - 1, particle.lon]) ** 2
            vel_horizontal = (u_horizontal + v_horizontal) ** (1/2)
            #
            u_star = (vel_horizontal * k) / ((math.log(fieldset.e3t[time, fieldset.mbathy[time, particle.depth, particle.lat, particle.lon] -1, particle.lat, particle.lon] / 2) / z_star))
            # Here tau is the bottom friction parameter estimated from (u_starr)^2 x density
            tau = ((u_star) ** 2) * 1024
            #
            #
            if tau >= threshold: # for colloids and marine particles
                frac_value = ParcelsRandom.randint(0,10)
                if frac_value >= 3:
                    particle.status = 2
                else:
                    particle.status = 3    
            #    
            else:  # for particles staying at the bottom
                particle.status = 4
    #
    #### OTHERS ####
    #def export(particle,fieldset,time):
    #    if particle.lat<48.7 and particle.lon < -124.66:
    #        particle.status = 7
    #
    def CheckOutOfBounds(particle, fieldset, time):
        if particle.state == StatusCode.ErrorOutOfBounds:    
            particle.delete()
    #        
    def KeepInOcean(particle, fieldset, time):
        if particle.state == StatusCode.ErrorThroughSurface:
            particle.depth = 0.0
            particle.state = StatusCode.Success             
    ##########################################################
    #
    #
    pset_states.execute([PBDEs_states, PBDEs_forms, Advection, turb_mix, Displacement, resuspension, CheckOutOfBounds, export, KeepInOcean],
                runtime=duration, 
                dt=dt_resolution,# dt resolution in seconds
                output_file=pset_states.ParticleFile(name=outfile_states, outputdt=timedelta(hours=odt)))
                #
    ds = xr.open_zarr(outfile_states)
    depth = ds.z/ds.fact           
    #
    return outfile_states