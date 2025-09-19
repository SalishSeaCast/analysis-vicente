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

sys.path.append('/home/vicentev/projects/def-allen/vicentev/analysis-vicente/OP_nibi')
from OP_functions_nibi import *


def timings(start_time, sim_length, number_outputs):
    month_days = 30 # 30  # number of days to release particles
    data_length = max(sim_length, 1)
    duration = datetime.timedelta(days=sim_length)
    delta_t = 5  # seconds
    release_particles_every = 900  # seconds

    number_particles = int(min(sim_length, month_days) * 86400 / release_particles_every)
    print("number_particles", number_particles)

    output_interval = datetime.timedelta(seconds=sim_length * 86400 / number_outputs)
    print('output_interval', output_interval)

    return start_time, data_length, duration, delta_t, release_particles_every, number_particles, output_interval


def name_outfile(year, month, sim_length, string):
    path = '/home/vicentev/scratch/vicentev/Simulations_Runs/'
    fn = f'PBDEs_01{month}{year}_run_{sim_length}_days_' + string + '.zarr'
    return os.path.join(path, fn)


def set_fieldsets_and_constants(start_time, data_length, delta_t):
    constants = {}
    constants['Iona_clat'] = [49.195045]
    constants['Iona_clon'] = [-123.301956]
    constants['Iona_z'] = 70
    constants['fraction_colloidal'] = 0.25

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
    field_set.add_constant('sinkvel_sewage', 12.84 * dt_h) # 12.84 m / hr --> 0.0035 m/s
    field_set.add_constant('sinkvel_marine', 5.52 * dt_h) # 2 m / hr    # 5.52 m / hr   # 12 m/hr
    ratio_MC = 0.2 # 0.08 #0.065 #0.1 #0.2 #0.4 # 0.012 # Ratio between Dissolved and Particulate PBDEs in the water column (Based on Sun et al., 2023)
    abso = (0.052 / 86400) #(0.024 / 86400) ##(0.038 / 86400)  # Colloidal/Dissolved → Attached to Marine Particle /s
    deso_s = (abso / ratio_MC) # Sewage Particle → Colloidal/Dissolved /s
    deso_m = (abso / ratio_MC) # Marine Particle → Colloidal/Dissolved /s
    deso_sed = deso_m
    abso_sed = deso_sed * 30. / 70

    field_set.add_constant('abso_probability', 1 - np.exp(-abso * delta_t))
    field_set.add_constant('deso_s_probability', 1 - np.exp(-deso_s * delta_t))
    field_set.add_constant('deso_m_probability', 1 - np.exp(-deso_m * delta_t))
    field_set.add_constant('deso_sed_probability', 1 - np.exp(-deso_sed * delta_t))
    field_set.add_constant('abso_sed_probability', 1 - np.exp(-abso_sed * delta_t))

    deg2met = 111319.5
    latT = 0.6495
    field_set.add_constant('u_deg2mps', deg2met * latT)
    field_set.add_constant('v_deg2mps', deg2met)

    kappa = 0.42
    zo, rho = 0.07, 1024
    field_set.add_constant('log_z_star', np.log(zo))
    cdmin, cdmax = 0.0075, 2
    field_set.add_constant('lowere3t_o2', zo * np.exp(kappa / np.sqrt(cdmax)))
    field_set.add_constant('uppere3t_o2', zo * np.exp(kappa / np.sqrt(cdmin)))

    tau_crit = 0.01#0.02 #0.025 # 0.01 #0.05 # 0.25
    field_set.add_constant('tau_constant', tau_crit / ((kappa ** 2) * rho))
    field_set.add_constant('tau_constant_lower', tau_crit / (rho * cdmax))
    field_set.add_constant('tau_constant_upper', tau_crit / (rho * cdmin))
    #
    field_set.add_constant('dx_lat', 0.00189/2)
    field_set.add_constant('dx_lon', 0.00519/2)
    field_set.add_constant('dy_lat', 0.00393/2)
    field_set.add_constant('dy_lon', -0.00334/2)

    return field_set, constants


def PBDEs_OP_run(year, month, day, sim_length, number_outputs, string,
                 restart=False, restart_filename=None):

    if restart and restart_filename is not None:
        ds = xr.open_zarr(restart_filename)
        last_time = np.nanmax(ds.time.values)
        start_time = pd.to_datetime(last_time).to_pydatetime()
        print(f"Restarting from {start_time} based on {restart_filename}")
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
                                    initial=np.arange(0, release_particles_every * number_particles, release_particles_every))
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
        restart_output_dir = "/home/vicentev/scratch/vicentev/Simulations_Runs/RESTART_Runs"
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
    print("RUNNING PBDEs_OP_run! :D")

    if len(sys.argv) == 5:
        sim_length = int(sys.argv[1])
        number_outputs = int(sys.argv[2])
        name_extension = str(sys.argv[3])
        restart_basename = str(sys.argv[4])
        path = '/home/vicentev/scratch/vicentev/Simulations_Runs'
        restart_filename = os.path.join(path, restart_basename)
        PBDEs_OP_run(None, None, None, sim_length, number_outputs,
                     name_extension, restart=True, restart_filename=restart_filename)

    elif len(sys.argv) == 7:
        year = int(sys.argv[1])
        month = int(sys.argv[2])
        day = int(sys.argv[3])
        sim_length = int(sys.argv[4])
        number_outputs = int(sys.argv[5])
        name_extension = str(sys.argv[6])
        PBDEs_OP_run(year, month, day, sim_length, number_outputs, name_extension)

    else:
        print("Invalid number of arguments!")
        print("Usage (restart): python -m Model_Driver_RESTART sim_length number_outputs name_extension restart_file.zarr")
        print("Usage (normal): python -m Model_Driver_RESTART year month day sim_length number_outputs name_extension")


    #
    ## How to run in the terminal as normal:
    #
    # python -m Model_Driver start_year start_month start_day length_sim_in_days number_outputs name_extension
    #
    ## How to run with restart file:
    #
    # python -m Model_Driver length_sim_in_days number_outputs name_extension_new restart_file