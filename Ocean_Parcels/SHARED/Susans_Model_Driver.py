import calendar
import datetime
import importlib 
import numpy as np
import os
import sys
import xarray as xr

from parcels import Field, FieldSet, ParticleSet,Variable, JITParticle

import OP_functions_shared as OP
import PBDEs_OP_Kernels_shared as PBDE

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def timings(year, month, day, sim_length, number_outputs):
    start_time = datetime.datetime(year, month, day)
    month_days = 30 # number of days to release particles
    data_length = max(sim_length, 1)
    duration = datetime.timedelta(days=sim_length)
    delta_t = 5 # s
    release_particles_every = 900 # s

    number_particles = int(min(sim_length, month_days) * 86400 / release_particles_every)
    print (number_particles)

    output_interval = datetime.timedelta(seconds=sim_length * 86400 / number_outputs)
#    output_interval = datetime.timedelta(seconds=5)
    print ('output_interval', output_interval)

    return (start_time, data_length, duration, delta_t, release_particles_every, number_particles, output_interval)


def name_outfile(year, month, day, sim_length):
    path = '/home/sallen/MEOPAR/ANALYSIS/analysis-vicente/Ocean_Parcels/SHARED/'
    print (year, month, sim_length)
    fn = f'NewRunInfoCorBot_for_{day}-{month}-{year}_run_{sim_length}_days.zarr'
    return os.path.join(path, fn)


#### CREATING FIELDSETS and Setting Constants ####
def set_fieldsets_and_constants(start_time, data_length, delta_t):

    constants = {}

    # Iona Outfall Location
    constants['Iona_clat'] = [49.195045]
    constants['Iona_clon'] = [-123.301956]
    constants['Iona_z'] = 70 # m
    # Iona output sewage vs colloidal
    constants['fraction_colloidal'] = 0.25 
    
    # Velocities
    varlist = ['U', 'V', 'W']
    filenames, variables = OP.filename_set(start_time, data_length, varlist)
    dimensions = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw','time': 'time_counter'}
    field_set = FieldSet.from_nemo(filenames, variables, dimensions, allow_time_extrapolation=True, chunksize='auto',
                                  )
    
    # Vertical Variables and Depth Related
    varlist=['Kz', 'totdepth', 'e3t', 'ssh', 'tmask', 'umask', 'vmask', 'fmask']
    filenames, variables = OP.filename_set(start_time, data_length, varlist)
    
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
    e3t = Field.from_netcdf(filenames['e3t'], variables['e3t'], dimensions, allow_time_extrapolation=True, chunksize='auto')
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
    field_set.add_constant('sinkvel_sewage', 500/86400.) # m/hr * dt --> to seconds --> ~ 500 m/d 
    field_set.add_constant('sinkvel_marine', 250/86400.) # m/hr * dt --> to seconds --> ~ 250 m/d 

    abso = 0.1 / 86400 # Colloidal/Dissolved → Attached to Marine Particle /s  (1/10 days)
    deso_s = 1.6 / 86400 # Sewage Particle → Colloidal/Dissolved /s
    deso_m = abso / 0.2 # Marine Particle → Colloidal/Dissolved /s ( 1 day)
    deso_sed = 1 / (7. * 86400) # slower than watercolumn (7 days)
    abso_sed = deso_sed * 30. / 70 # in the sediments, easier to find marine particles to bind to, 30/70 is ratio of suspended materials (1/14.6 days)
    sediment_burying = 1. / (10000 * 365 * 86400) # Particles get buried by sediment
    field_set.add_constant('abso_probability', 1 - np.exp(-abso * delta_t))
    field_set.add_constant('deso_s_probability', 1 - np.exp(-deso_s * delta_t))
    field_set.add_constant('deso_m_probability', 1 - np.exp(-deso_m * delta_t))
    field_set.add_constant('deso_sed_probability', 1 - np.exp(-deso_sed * delta_t))
    field_set.add_constant('abso_sed_probability', 1 - np.exp(-abso_sed * delta_t))
    field_set.add_constant('sediment_burying_probability', 1 - np.exp(-sediment_burying * delta_t))
    print (field_set.abso_probability, field_set.deso_s_probability, field_set.deso_m_probability)
    print (field_set.abso_sed_probability, field_set.deso_sed_probability, field_set.sediment_burying_probability)
    
    # conversion factors
    deg2met = 111319.5
    latT = 0.6495
    field_set.add_constant('u_deg2mps', deg2met*latT)
    field_set.add_constant('v_deg2mps', deg2met)
    
    kappa = 0.42
    zo, rho = 0.07, 1024                              # from SalishSeaCast
    field_set.add_constant('log_z_star', np.log(zo))
    cdmin, cdmax = 0.0075, 2                          # from SalishSeaCast
    field_set.add_constant('lowere3t_o2', zo * np.exp(kappa / np.sqrt(cdmax)))
    field_set.add_constant('uppere3t_o2', zo * np.exp(kappa / np.sqrt(cdmin)))

    tau_crit = 0.005 # halved
    tau_bury_crit = 0.8
    field_set.add_constant('tau_constant', tau_crit / ((kappa ** 2) * rho))
    field_set.add_constant('tau_constant_lower', tau_crit / (rho * cdmax))
    field_set.add_constant('tau_constant_upper', tau_crit / (rho * cdmin))
    field_set.add_constant('tau_bury_constant', tau_bury_crit / ((kappa ** 2) * rho))
    field_set.add_constant('tau_bury_constant_lower', tau_bury_crit / (rho * cdmax))
    field_set.add_constant('tau_bury_constant_upper', tau_bury_crit / (rho * cdmin))

    print (field_set.tau_constant, field_set.tau_bury_constant)
    print (field_set.tau_constant_lower, field_set.tau_constant_upper)

    field_set.add_constant('dx_lat', 0.00189/2)
    field_set.add_constant('dx_lon', 0.00519/2)
    field_set.add_constant('dy_lat', 0.00393/2)
    field_set.add_constant('dy_lon', -0.00334/2)
    
    return field_set, constants


def PBDEs_OP_run(year, month, day, sim_length, number_outputs, newstart=True, input_file=None):

    # Set-up Run
    (start_time, data_length, duration, delta_t, 
         release_particles_every, number_particles, output_interval) = timings(year, month, day, sim_length, number_outputs)

    field_set, constants = set_fieldsets_and_constants(start_time, data_length, delta_t)

    outfile_states = name_outfile(year, month, day, sim_length)

    # Set-up Ocean Parcels
    class MPParticle(JITParticle):
        status = Variable('status', initial=(np.random.rand(number_particles) >
                                             constants['fraction_colloidal']).astype(int) - 2)
        vvl_factor = Variable('fact', initial=1)
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
        release_time = Variable('release_time', 
                        initial=np.arange(0, release_particles_every*number_particles, release_particles_every))
    class MPParticle_R(JITParticle):
        release_time = Variable('release_time')
        status = Variable('status')
        fact = Variable('fact')

    print (newstart)
    if newstart:
        pset_states = ParticleSet(field_set, pclass=MPParticle, lon=constants['Iona_clon']*np.ones(number_particles), 
                          depth=constants['Iona_z']*np.ones(number_particles), 
                              lat = constants['Iona_clat']*np.ones(number_particles))
    else:
        input_file = input_file
        print (input_file)
        print (xr.open_dataset(input_file))
        pset_states = ParticleSet.from_particlefile(field_set, pclass=MPParticle_R, filename=input_file)

    output_file = pset_states.ParticleFile(name=outfile_states, outputdt=output_interval)
    
    KE = (pset_states.Kernel(PBDE.PBDEs_states) + pset_states.Kernel(PBDE.Sinking) 
      + pset_states.Kernel(PBDE.Advection)
      + pset_states.Kernel(PBDE.turb_mix) + pset_states.Kernel(PBDE.resuspension)
      + pset_states.Kernel(PBDE.CheckOutOfBounds) + pset_states.Kernel(PBDE.export)
      + pset_states.Kernel(PBDE.KeepInOcean)
     )

    # Run!
    pset_states.execute(KE, runtime=duration, dt=delta_t, output_file=output_file)
    

if __name__ == "__main__":
    print("RUNNING PBDEs_OP_run! :D")
    # Input from the terminal
    year = int(sys.argv[1])  # Example: 2022
    month = int(sys.argv[2]) # Integer between 1 and 12
    day = int(sys.argv[3]) # Interger between 1 and 31
    sim_length = int(sys.argv[4]) 
    number_outputs = int(sys.argv[5])
    if sys.argv[6] == 'True':
        newstart = True
        input_file = None
    else:
        newstart = False
        input_file = sys.argv[7]

    PBDEs_OP_run(year, month, day, sim_length, number_outputs, newstart, input_file)
    #
    ## How to run in the terminal:
    # python -m Susans_Model_Driver start_year start_month start_day length_sim_in_days number_outputs
