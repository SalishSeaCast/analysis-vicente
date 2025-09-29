import calendar
import datetime
import importlib 
import numpy as np
import os
import sys
from datetime import timedelta
import xarray as xr
import pandas as pd

from parcels import Field, FieldSet, ParticleSet, Variable, JITParticle

sys.path.append('/home/vicentev/projects/def-allen/vicentev/analysis-vicente/OP_nibi/base_OP')
from OP_functions_nibi import *

# Defining parameters from yaml file:
yaml_file_name = sys.argv[1]
config_yaml = [os.path.join('/home/vicentev/projects/def-allen/vicentev/analysis-vicente/OP_nibi/config_files/', yaml_file_name)]
param = load_config(config_yaml)
#Timing definitions
start_simulation = datetime.datetime(param['release_params']['start_sim_year'], param['release_params']['start_sim_month'], param['release_params']['start_sim_day']) #Start date simulation
start_release = datetime.datetime(param['release_params']['year_start_release'], param['release_params']['month_start_release'], param['release_params']['day_start_release']) #Start date release
day_release = param['release_params']['days_of_release'] # how many days to release particles
release_freq = param['release_params']['release_particles_freq'] # release frequency in seconds
length = param['release_params']['simulation_length'] # Simulation length in days
delta_t = param['release_params']['delta_t'] # Processes resolution in seconds
n_outputs = param['release_params']['number_outputs'] # number of output observations
#
# Iona Constants Definitions
lon_iona = param['constants']['lon_iona'] # Longitude coordinate
lat_iona = param['constants']['lat_iona'] # Latitude coordinate
depth_iona = param['constants']['depth_iona'] # Depth of release
fraction_colloidal = param['constants']['fraction_colloidal'] # fraction of particles released in colloidal phase
#
# Particles Features
vel_sewage = param['particles_features']['sinking_vel_sewage'] # sinking vel of sewage particles
vel_marine = param['particles_features']['sinking_vel_marine'] # sinking vel of marine particles
absorption = param['particles_features']['absorption'] # absorption of colloidal to marine particles
ratio_marine_colloidal = param['particles_features']['ratio_marine_colloidal'] # ratio between colloidal and marine particles in the WC
fraction_sediment = param['particles_features']['fraction_sediment'] # fraction of colloidal to marine particles in the sediment
#
# Grid Parameters
deg2met_value = param['grid_params']['deg2met'] # conversion from degrees to meters
latT_value = param['grid_params']['latT'] 
dx_lat_value = param['grid_params']['dx_lat']
dx_lon_value = param['grid_params']['dx_lon']
dy_lat_value = param['grid_params']['dy_lat']
dy_lon_value = param['grid_params']['dy_lon']
#
# Resuspension Parameters
kappa_value = param['resuspension_params']['kappa']
zo_value = param['resuspension_params']['zo']
rho_value = param['resuspension_params']['rho']
cdmin_value = param['resuspension_params']['cdmin']
cdmax_value = param['resuspension_params']['cdmax']
tau_crit_value = param['resuspension_params']['tau_critical'] # critical resuspension tau value
#
# Name Extension Simulation
extension = param['name_extension']
#
##### Initilization time for release in seconds
seconds_initial = (start_release - start_simulation).total_seconds()
#####

def timings(start_time, sim_length, number_outputs):
    data_length = max(sim_length, 1)
    duration = datetime.timedelta(days=sim_length)

    number_particles = int(min(sim_length, day_release) * 86400 / release_freq)
    print("number_particles", number_particles)

    output_interval = datetime.timedelta(seconds=sim_length * 86400 / number_outputs)
    print('output_interval', output_interval)

    return start_time, data_length, duration, delta_t, release_freq, number_particles, output_interval


def name_outfile(year, month , sim_length, string):
    path = param['simulations_output_dir']
    fn = f'PBDEs_01{month}{year}_run_{sim_length}_days_' + string + '.zarr'
    return os.path.join(path, fn)


def set_fieldsets_and_constants(start_time, data_length, delta_t):
    constants = {}
    constants['Iona_clat'] = [lat_iona]
    constants['Iona_clon'] = [lon_iona]
    constants['Iona_z'] = [depth_iona]
    constants['fraction_colloidal'] = fraction_colloidal

    varlist = ['U', 'V', 'W']
    filenames, variables = filename_set(start_time, data_length, varlist)
    dimensions = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}
    field_set = FieldSet.from_nemo(filenames, variables, dimensions, allow_time_extrapolation=True, chunksize='auto')
    #
    varlist=['Kz', 'totdepth', 'cell_size', 'ssh', 'tmask', 'umask', 'vmask', 'fmask']
    filenames, variables = filename_set(start_time, data_length, varlist)
    #
    # 2-D, no time
    dimensions = {'lon': 'glamt', 'lat': 'gphit'}
    TD = Field.from_netcdf(filenames['totdepth'], variables['totdepth'], dimensions, allow_time_extrapolation=True, chunksize='auto')
    field_set.add_field(TD)
    # 2-D, with time
    dimensions = {'lon': 'glamt', 'lat': 'gphit', 'time': 'time_counter'}
    SSH = Field.from_netcdf(filenames['ssh'], variables['ssh'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    field_set.add_field(SSH)
    # 3-D on W
    dimensions = {'lon': 'glamt', 'lat': 'gphit', 'depth': 'depthw','time': 'time_counter'}
    Kz = Field.from_netcdf(filenames['Kz'], variables['Kz'], dimensions, allow_time_extrapolation=True, chunksize='auto')
    field_set.add_field(Kz)
    # 3-D on T
    dimensions = {'lon': 'glamt', 'lat': 'gphit', 'depth': 'deptht','time': 'time_counter'}
    e3t = Field.from_netcdf(filenames['cell_size'], variables['cell_size'], dimensions, allow_time_extrapolation=True, chunksize='auto')
    field_set.add_field(e3t)

    dimensions = {'lon': 'glamt', 'lat': 'gphit', 'depth': 'deptht','time': 't'}
    tmask = Field.from_netcdf(filenames['tmask'], variables['tmask'], dimensions, chunksize='auto')
    field_set.add_field(tmask)
    
    dimensions = {'lon': 'glamu', 'lat': 'gphiu', 'depth': 'deptht','time': 't'}
    umask = Field.from_netcdf(filenames['umask'], variables['umask'], dimensions, chunksize='auto')
    field_set.add_field(umask)
    
    dimensions = {'lon': 'glamv', 'lat': 'gphiv', 'depth': 'deptht','time': 't'}
    vmask = Field.from_netcdf(filenames['vmask'], variables['vmask'], dimensions, chunksize='auto')
    field_set.add_field(vmask)
    
    dimensions = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'deptht','time': 't'}
    fmask = Field.from_netcdf(filenames['fmask'], variables['fmask'], dimensions, chunksize='auto')
    field_set.add_field(fmask)

    dt_h = 1 / 3600.
    frac_sed = fraction_sediment #30. / 70
    field_set.add_constant('sinkvel_sewage', vel_sewage * dt_h) # 12.84 m / hr --> 0.0035 m/s
    field_set.add_constant('sinkvel_marine', vel_marine * dt_h) # 2 m / hr    # 5.52 m / hr   # 12 m/hr
    ratio_MC = ratio_marine_colloidal # 0.08 #0.065 #0.1 #0.2 #0.4 # 0.012 # Ratio between Dissolved and Particulate PBDEs in the water column (Based on Sun et al., 2023)
    abso = (absorption / 86400) #(0.024 / 86400) ##(0.038 / 86400)  # Colloidal/Dissolved → Attached to Marine Particle /s
    deso_s = (abso / ratio_MC) # Sewage Particle → Colloidal/Dissolved /s
    deso_m = (abso / ratio_MC) # Marine Particle → Colloidal/Dissolved /s
    deso_sed = deso_m
    abso_sed = deso_sed * frac_sed

    field_set.add_constant('abso_probability', 1 - np.exp(-abso * delta_t))
    field_set.add_constant('deso_s_probability', 1 - np.exp(-deso_s * delta_t))
    field_set.add_constant('deso_m_probability', 1 - np.exp(-deso_m * delta_t))
    field_set.add_constant('deso_sed_probability', 1 - np.exp(-deso_sed * delta_t))
    field_set.add_constant('abso_sed_probability', 1 - np.exp(-abso_sed * delta_t))

    deg2met = deg2met_value
    latT = latT_value
    field_set.add_constant('u_deg2mps', deg2met * latT)
    field_set.add_constant('v_deg2mps', deg2met)

    kappa = kappa_value
    zo, rho = zo_value, rho_value
    field_set.add_constant('log_z_star', np.log(zo))
    cdmin, cdmax = cdmin_value, cdmax_value
    field_set.add_constant('lowere3t_o2', zo * np.exp(kappa / np.sqrt(cdmax)))
    field_set.add_constant('uppere3t_o2', zo * np.exp(kappa / np.sqrt(cdmin)))

    tau_crit = tau_crit_value#0.01#0.02 #0.025 # 0.01 #0.05 # 0.25
    field_set.add_constant('tau_constant', tau_crit / ((kappa ** 2) * rho))
    field_set.add_constant('tau_constant_lower', tau_crit / (rho * cdmax))
    field_set.add_constant('tau_constant_upper', tau_crit / (rho * cdmin))
    #
    field_set.add_constant('dx_lat', dx_lat_value)
    field_set.add_constant('dx_lon', dx_lon_value)
    field_set.add_constant('dy_lat', dy_lat_value)
    field_set.add_constant('dy_lon', dy_lon_value)

    return field_set, constants


def PBDEs_OP_run(year , month, day, sim_length, number_outputs, string,
                 restart=False, restart_filename=None):

    if restart and restart_filename is not None:
        ds = xr.open_zarr(restart_filename)
        last_time = np.nanmax(ds.time.values)
        start_time = pd.to_datetime(last_time).to_pydatetime()
        ds.close()
    else:
        start_time = datetime.datetime(year, month, day)

    (start_time, data_length, duration, delta_t, release_particles_every,
     number_particles, output_interval) = timings(start_time, sim_length, number_outputs)

    field_set, constants = set_fieldsets_and_constants(start_time, data_length, delta_t)

    if restart and restart_filename is not None:
        class MPParticle(JITParticle):
            status = Variable('status')
            fact = Variable('fact')
            release_time = Variable('release_time')
            #
            tmask = Variable('tmask')
            umask = Variable('umask')
            vmask = Variable('vmask')
            fmask = Variable('fmask')
            stuck = Variable('stuck')
            H_vel_2 = Variable('H_vel_2')
            crit = Variable('crit')
            bat_particle = Variable('bat_particle')
            uvalue = Variable('uvalue')
            vvalue = Variable('vvalue')
            wvalue = Variable('wvalue')
            e3t = Variable('e3t')
            totaldepth = Variable('totaldepth')
    else:
        class MPParticle(JITParticle):
            status = Variable('status', initial=(np.random.rand(number_particles) >
                                                constants['fraction_colloidal']).astype(int) - 2)
            fact = Variable('fact', initial=1)
            release_time = Variable('release_time',
                                    initial=np.arange(seconds_initial, seconds_initial + (release_freq * number_particles), release_freq))
            #
            tmask = Variable('tmask', initial=1)
            umask = Variable('umask', initial=1)
            vmask = Variable('vmask', initial=1)
            fmask = Variable('fmask', initial=1)
            stuck = Variable('stuck', initial=0)
            H_vel_2 = Variable('H_vel_2', initial=0)
            crit = Variable('crit', initial=0)
            bat_particle = Variable('bat_particle', initial=0)
            uvalue = Variable('uvalue', initial=0)
            vvalue = Variable('vvalue', initial=0)
            wvalue = Variable('wvalue', initial=0)
            e3t = Variable('e3t', initial=0)
            totaldepth = Variable('totaldepth', initial=0)          


    if restart and restart_filename is not None:
        pset_states = ParticleSet.from_particlefile(fieldset=field_set, pclass=MPParticle, filename=restart_filename)
        restart_output_dir = param['restart_output_dir']
        restart_basename = os.path.basename(restart_filename).replace('.zarr', f'_restart_{sim_length}_days_{name_extension}.zarr')
        outfile_states = os.path.join(restart_output_dir, restart_basename)
    else:
        pset_states = ParticleSet(field_set, pclass=MPParticle,
                                  lon=constants['Iona_clon'] * np.ones(number_particles),
                                  depth=constants['Iona_z'] * np.ones(number_particles),
                                  lat=constants['Iona_clat'] * np.ones(number_particles))
        outfile_states = name_outfile(year, month, sim_length, string)

    output_file = pset_states.ParticleFile(name=outfile_states, outputdt=output_interval)

    KE = (pset_states.Kernel(PBDEs_states) + pset_states.Kernel(Sinking) +
          pset_states.Kernel(Advection) + pset_states.Kernel(turb_mix) +
          pset_states.Kernel(resuspension)  + pset_states.Kernel(export_JdF) + 
          pset_states.Kernel(export_Js) + pset_states.Kernel(CheckOutOfBounds) + pset_states.Kernel(KeepInOcean))

    pset_states.execute(KE, runtime=duration, dt=delta_t, output_file=output_file)


if __name__ == "__main__":
    #
    if param['RESTART_true_or_false'] == True:
        sim_length = param['release_params']['simulation_length']
        number_outputs = param['release_params']['number_outputs']
        name_extension = param['name_extension']
        restart_basename = param['restart_filename']
        #path = '/home/vicentev/scratch/vicentev/Simulations_Runs'
        #restart_filename = os.path.join(path, restart_basename)
        restart_filename = param['restart_filename']
        print('Running Simulation from RESTART file ', restart_filename)
        PBDEs_OP_run(None, None, None, sim_length, number_outputs,
                     name_extension, restart=True, restart_filename=restart_filename)

    elif param['RESTART_true_or_false'] == False:
        year = param['release_params']['start_sim_year']
        month = param['release_params']['start_sim_month']
        day = param['release_params']['start_sim_day']
        sim_length = param['release_params']['simulation_length']
        number_outputs = param['release_params']['number_outputs']
        name_extension = param['name_extension']
        print('Running new Simulation from yaml file')
        PBDEs_OP_run(year, month, day, sim_length, number_outputs, name_extension)
        #  The yaml file should include this variables to initialize the simulations. 
        # Add some others for other purposes; times for release ...
    else:
        print("Something went wrong! Check yaml file! :O")


    #
    ## How to run in the terminal:
    #
    # python -m Model_Driver yaml_file.yaml
    #
    # 1) Do yaml files to input parameters and run smoothly from the terminal and set output.txt file (DONE)
    # 2) Do test runs to check later parallel running