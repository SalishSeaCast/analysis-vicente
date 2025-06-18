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

sys.path.append('/ocean/vvalenzuela/MOAD/Ocean_Parcels')
from OP_functions import *


def timings(start_time, sim_length, number_outputs):
    month_days = 30  # number of days to release particles
    data_length = max(sim_length, 1)
    duration = datetime.timedelta(days=sim_length)
    delta_t = 5  # seconds
    release_particles_every = 1800  # seconds

    number_particles = int(min(sim_length, month_days) * 86400 / release_particles_every)
    print("number_particles", number_particles)

    output_interval = datetime.timedelta(seconds=sim_length * 86400 / number_outputs)
    print('output_interval', output_interval)

    return start_time, data_length, duration, delta_t, release_particles_every, number_particles, output_interval


def name_outfile(year, month, sim_length, string):
    path = '/home/vvalenzuela/MOAD/Ocean_Parcels/results/Simulations_runs'
    fn = f'PBDE_particles_for_01{month}{year}_run_{sim_length}_days_' + string + '.zarr'
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

    varlist = ['Kz', 'totdepth', 'cell_size', 'ssh']
    filenames, variables = filename_set(start_time, data_length, varlist)

    field_set.add_field(Field.from_netcdf(filenames['totdepth'], variables['totdepth'],
                                          {'lon': 'glamt', 'lat': 'gphit'},
                                          allow_time_extrapolation=True, chunksize='auto'))
    field_set.add_field(Field.from_netcdf(filenames['ssh'], variables['ssh'],
                                          {'lon': 'glamt', 'lat': 'gphit', 'time': 'time_counter'},
                                          allow_time_extrapolation=True, chunksize='auto'))
    field_set.add_field(Field.from_netcdf(filenames['Kz'], variables['Kz'],
                                          {'lon': 'glamt', 'lat': 'gphit', 'depth': 'depthw', 'time': 'time_counter'},
                                          allow_time_extrapolation=True, chunksize='auto'))
    field_set.add_field(Field.from_netcdf(filenames['cell_size'], variables['cell_size'],
                                          {'lon': 'glamt', 'lat': 'gphit', 'depth': 'deptht', 'time': 'time_counter'},
                                          allow_time_extrapolation=True, chunksize='auto'))

    dt_h = 1 / 3600.
    field_set.add_constant('sinkvel_sewage', 12.84 * dt_h) # 12.84 m / hr
    field_set.add_constant('sinkvel_marine', 5.52 * dt_h) # 2 m / hr    # 5.52 m / hr   # 9 m/hr
    ratio_MP = 0.13
    abso = (0.038 / 86400)  # Colloidal/Dissolved → Attached to Marine Particle /s
    deso_s = (abso / ratio_MP) # Sewage Particle → Colloidal/Dissolved /s
    deso_m = (abso / ratio_MP) # Marine Particle → Colloidal/Dissolved /s
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

    tau_crit = 0.01 #0.0025 #0.04
    field_set.add_constant('tau_constant', tau_crit / ((kappa ** 2) * rho))
    field_set.add_constant('tau_constant_lower', tau_crit / (rho * cdmax))
    field_set.add_constant('tau_constant_upper', tau_crit / (rho * cdmin))

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
    else:
        class MPParticle(JITParticle):
            status = Variable('status', initial=(np.random.rand(number_particles) >
                                                constants['fraction_colloidal']).astype(int) - 2)
            fact = Variable('fact', initial=1)
            release_time = Variable('release_time',
                                    initial=np.arange(0, release_particles_every * number_particles, release_particles_every))


    if restart and restart_filename is not None:
        pset_states = ParticleSet.from_particlefile(fieldset=field_set, pclass=MPParticle, filename=restart_filename)
        restart_output_dir = "/home/vvalenzuela/MOAD/Ocean_Parcels/results/Simulations_runs/RESTART_runs"
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
          pset_states.Kernel(resuspension) + pset_states.Kernel(CheckOutOfBounds) +
          pset_states.Kernel(export) + pset_states.Kernel(KeepInOcean))

    pset_states.execute(KE, runtime=duration, dt=delta_t, output_file=output_file)


if __name__ == "__main__":
    print("RUNNING PBDEs_OP_run! :D")

    if len(sys.argv) == 5:
        sim_length = int(sys.argv[1])
        number_outputs = int(sys.argv[2])
        name_extension = str(sys.argv[3])
        restart_basename = str(sys.argv[4])
        path = '/home/vvalenzuela/MOAD/Ocean_Parcels/results/Simulations_runs'
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