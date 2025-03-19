#
### Ocean Parcels Kernels for PBDEs Simulations ###    
#
#
def PBDEs_states(particle, fieldset, time):    
    #
    if particle.release_time <= 0:
        particle.release_time = particle.time
    # Particle.release_time is the internal time where particles are being added     
    # 
    # the next if statement change to status 1 or 2 at the time of release without overwritting older particles
    if particle.time == particle.release_time:
        if ParcelsRandom.uniform(0, 1) < 0.75:  # 75% chance for status 1
            particle.status = 1  # Sewage Particle
        else:
            particle.status = 2  # Colloidal/Dissolved PBDEs
    else:
    #
        random_value = ParcelsRandom.random()
        # Absorption and desoprtion (per hour)
        abso = 0.7 / 24   # Colloidal/Dissolved → Attached to Marine Particle
        deso_s = 3.2 / 24 # Sewage Particle → Colloidal/Dissolved
        deso_m = 1.8 / 24 # Marine Particle → Colloidal/Dissolved

        # Status updates
        if particle.status == 2 and random_value < 1 - math.exp(-abso * particle.dt_h):
            particle.status = 3  # Becomes attached to Marine Particle
        elif particle.status == 1 and random_value < 1 - math.exp(-deso_s * particle.dt_h):
            particle.status = 2  # Becomes Colloidal/Dissolved
        elif particle.status == 3 and random_value < 1 - math.exp(-deso_m * particle.dt_h):
            particle.status = 2  # Returns to Colloidal/Dissolved
    #
#### PBDEs states sinking velocities features ####
def PBDEs_forms(particle, fieldset, time):
    #print('PBDEs_forms kernel is running')
    #
    #dt_h = 1 / 3600
    if particle.status == 1:
        sinkvel = 50*(particle.dt_h) # m/hr * dt --> to seconds
        particle.depth += sinkvel * particle.dt
    #Sewage Particles sink fast        
    elif particle.status == 2:
        sinkvel = 0.0
        particle.depth += sinkvel * particle.dt
    # Colloids just float around and move with advection
    elif particle.status == 3:
        sinkvel = 10*(particle.dt_h) # m/hr * dt --> to seconds
        particle.depth += sinkvel * particle.dt 
#
#### ADVECTION ####
def Advection(particle, fieldset, time):
    #print('Advection kernel is running') 
    # Advection for all PBDEs in status 1, 2 and 3
    if particle.status == 1 or particle.status == 2 or particle.status == 3: 
        ssh = fieldset.sossheig[time, particle.depth, particle.lat, particle.lon] #SSH(t) sea surface height
        sshn = fieldset.sossheig[time+particle.dt, particle.depth, particle.lat, particle.lon] #SSH(t+dt) sea surface height in the next time step
        td = fieldset.totaldepth[time, particle.depth, particle.lat, particle.lon]#Total_depth 
        particle.fact = (1+ssh/td)
        VVL = (sshn-ssh)*particle.depth/(td)
        #VVL = (sshn-ssh)*particle.depth/(td+ssh)
        #
        # calculate once and reuse
        dt_factor = 0.5 * particle.dt / particle.fact
        (u1, v1, w1) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
        lon1 = particle.lon + u1*.5*particle.dt
        lat1 = particle.lat + v1*.5*particle.dt
        dep1 = particle.depth + w1 * dt_factor
        #
        (u2, v2, w2) = fieldset.UVW[time + .5 * particle.dt, dep1, lat1, lon1]
        lon2 = particle.lon + u2*.5*particle.dt
        lat2 = particle.lat + v2*.5*particle.dt
        dep2 = particle.depth + w2 * dt_factor
        #
        (u3, v3, w3) = fieldset.UVW[time + .5 * particle.dt, dep2, lat2, lon2]
        lon3 = particle.lon + u3*particle.dt
        lat3 = particle.lat + v3*particle.dt
        dep3 = particle.depth + w3 * dt_factor
        #
        (u4, v4, w4) = fieldset.UVW[time + particle.dt, dep3, lat3, lon3]
        #
        wa = (w1 + 2*w2 + 2*w3 + w4) /6.
        particle.wa = wa* particle.dt
        particle_dlon = (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle_dlat = (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        particle_ddepth = particle.wa/particle.fact + VVL
        #
        if particle_ddepth + particle.depth < 0:
            particle_ddepth = - (2*particle.depth+particle_ddepth)
#
#### TURBULENT MIX ####
def turb_mix(particle,fieldset,time):
    #print('Turb_mix kernel is running')
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
    #print('Displacement kernel is running')
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
    # Only proceed if particle is at the seabed (status 4)
    if particle.status == 4:
        # Update e3t value only if the particle moved
        if (particle.depth != particle.prev_depth or 
            particle.lat != particle.prev_lat or 
            particle.lon != particle.prev_lon):
            #
            particle.e3t_val = fieldset.e3t[0, particle.depth, particle.lat, particle.lon]  # Update e3t
            particle.prev_depth = particle.depth
            particle.prev_lat = particle.lat
            particle.prev_lon = particle.lon  
        #
        bat_particle = particle.depth - 0.5 * particle.e3t_val  
        #
        # horizontal velocities in m/s  
        u_vel = fieldset.U[time, bat_particle, particle.lat, particle.lon] * particle.deg2met * particle.latT
        v_vel = fieldset.V[time, bat_particle, particle.lat, particle.lon] * particle.deg2met
        # squared horizontal velocity
        H_vel_2 = u_vel**2 + v_vel**2  
        particle.h_vel = H_vel_2  

        # Apply resuspension only if velocity exceeds threshold
        if H_vel_2 > 0.00029:
            log_e3t = math.log(particle.e3t_val * 0.5)  
            TAU = H_vel_2 * particle.k_constant * ((log_e3t - particle.log_z_star) ** (-2)) * 1024  
            particle.tau_values = TAU  

            # Resuspension probability (30% for marine particles, 70% for colloidal)
            if ParcelsRandom.uniform(0, 1) >= 0.3:
                particle.status = 2  # Colloidal/Dissolved PBDEs
            else:
                particle.status = 3  # Marine Particles
        else:
            particle.status = 4  # Particle remains at the seabed         
#
#### OTHERS ####
def export(particle,fieldset,time):
    if particle.lat<48.7 and particle.lon < -124.66:
        particle.status = 7

def CheckOutOfBounds(particle, fieldset, time):
    if particle.state == StatusCode.ErrorOutOfBounds:    
        particle.delete()
#        
def KeepInOcean(particle, fieldset, time):
    if particle.state == StatusCode.ErrorThroughSurface:
        particle.depth = 0.0
        particle.state = StatusCode.Success  