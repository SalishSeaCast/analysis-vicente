#
### Ocean Parcels Kernels for PBDEs Simulations ###    
#
# Force kernels to not move particles when not released!! <---- Not to move or change state
#
def PBDEs_states(particle, fieldset, time):    
    #   
    # Particle.release_time is the internal time where particles are being added (example: 0 3600 7200, etc)    
    #
    # INITIAL STATUS FORCED FROM OUTSIDE THE KERNEL
        #particle.status = particle.status_initial # initialized with 75% for status 1 and 25% for status 2 
    #
    #    
#    particle.age_h = (particle.time - particle.release_time) / 3600.0
    #    particle.release_time = particle.time
    # If the particle has just been released (for example, within one time step)
    # then ensure that its status reflects the release status.
#    if particle.age_h < particle.dt_h:
        # Optionally, you can force the status to equal release_status.
        # This way, if you're outputting at the release hour, you'll see the original value.
#        particle.status = particle.release_status
    #if particle.time < 10000:
        #particle.status = 0
    #
    #else:
    #    
    if (time > particle.release_time):
        if particle.status < 0:
            particle.status = - particle.status
        elif (particle.status < 4):
            random_value = ParcelsRandom.random()
        # Absorption and desoprtion (per hour)
        # Status updates
            if particle.status == 1 and random_value < fieldset.deso_s_probability:
                particle.status = 2  # Becomes Colloidal/Dissolved
            elif particle.status == 2 and random_value < fieldset.abso_probability:
                particle.status = 3  # Becomes attached to Marine Particle
            elif particle.status == 3 and random_value < fieldset.deso_m_probability:
                particle.status = 2  # Returns to Colloidal/Dissolved         
#### PBDEs states sinking velocities features ####
def PBDEs_forms(particle, fieldset, time):
    #print('PBDEs_forms kernel is running')
    if particle.status == 1:
        #particle_ddepth += fieldset.sinkvel_sewage * particle.dt
        particle.depth += fieldset.sinkvel_sewage * particle.dt
    #Sewage Particles sink fast        
    elif particle.status == 3:
        #particle_ddepth += fieldset.sinkvel_marine * particle.dt
        particle.depth += fieldset.sinkvel_marine * particle.dt 
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
    if particle.status == 1 or particle.status == 2 or particle.status == 3:
        """Vertical mixing"""
        #Vertical mixing
        td = fieldset.totaldepth[time, particle.depth, particle.lat, particle.lon] #Total_depth (can't share, unless I set it for the particle)
        if particle.depth + 0.5 / particle.fact > td: #Only calculate gradient of diffusion for particles deeper than 0.5 otherwise OP will check for particles outside the domain and remove it.
            Kzdz = 2 * (fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon] 
                        - fieldset.vert_eddy_diff[time, particle.depth-0.5/particle.fact, particle.lat, particle.lon]) #backwards difference 
        else: 
            Kzdz = 2 * (fieldset.vert_eddy_diff[time, particle.depth+0.5/particle.fact, particle.lat, particle.lon]
                        - fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon]) #forward difference 
        dgrad = Kzdz * particle.dt / particle.fact
        if particle.depth + (0.5 * dgrad) > 0 and particle.depth + (0.5 * dgrad) < td:
            Kz = fieldset.vert_eddy_diff[time, particle.depth+0.5*dgrad, particle.lat, particle.lon] #Vertical diffusivity SSC  
        else:
            Kz = 0 
        #
        Rr = ParcelsRandom.uniform(-1, 1)
        d_random = sqrt(3 * 2 * Kz * particle.dt) * Rr / particle.fact
        dzs = (dgrad + d_random)
        
        #Apply turbulent mixing.
        if dzs + particle_ddepth + particle.depth > td:
            particle.depth  = td # Get particles attached to the bottom when they reach it
            particle.status = 4
 #ADD LATER!!!           particle.e3t_val = fieldset.e3t[0, particle.depth, particle.lat+particle_dlat, particle.lon+particle_dlon]  # Update e3t
            #
        elif dzs + particle.depth + particle_ddepth < 0:
            particle_ddepth = -(dzs + particle.depth+particle_ddepth) #reflection on surface
        #
        else:
            particle_ddepth += dzs #apply mixing
#
#### VERTICAL DISPLACEMENT ####
#def Displacement(particle,fieldset,time):
    #print('Displacement kernel is running')
#    if particle.status == 1 or particle.status == 2 or particle.status == 3:
        #Apply turbulent mixing.
#        if dzs + particle_ddepth + particle.depth > td:
#            particle.depth  = td # Get particles attached to the bottom when they reach it
#            particle.status = 4
#            particle.e3t_val = fieldset.e3t[0, particle.depth, particle.lat+particle_dlat, particle.lon+particle_dlon]  # Update e3t
            #
#        elif dzs + particle.depth+ particle_ddepth < 0:
#            particle_ddepth = -(dzs + particle.depth+particle_ddepth) #reflection on surface
        #
#        else:
#            particle_ddepth += dzs #apply mixing
#
#### RESUSPENSION ####
# Try 2 min, 5 min, 10min for sensitivity :)
def resuspension(particle, fieldset, time):
    #particle.time_resuspension += math.fabs(particle.dt)
    #
    #if particle.time_resuspension >= particle.resuspension_interval:
        #
    #    particle.time_resuspension = 0 
        #
        # Only proceed if particle is at the seabed (status 4)
        if particle.status == 4:
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
            #if H_vel_2 > 0.00029:
                #log_e3t = math.log(particle.e3t_val * 0.5)  
            #TAU = H_vel_2 * particle.k_constant * ((particle.log_e3t - particle.log_z_star) ** (-2)) * 1024  
            #particle.tau_values = TAU
            # 
            # Intentar definir una aproximacion lineal para no tener que calcular el logaritmo!!! Te ahorrarias mucho tiempo de simulacion!!!
            factor = 1 + (fieldset.sossheig[time, particle.depth, particle.lat, particle.lon] / particle.depth) #SSH(t) sea surface height
            particle.log_e3t = math.log(particle.e3t_val * factor * 0.5)  
            # Taylor Expansion would need a for loop inside the Kernel :( ---> Heavier than the log function inside the Kernel!
            #
            # How to add the 10 minutes (for example) change in the probability of resuspending?
            if particle.tau_constant * (particle.log_e3t - particle.log_z_star) ** 2 >= H_vel_2: # new updated condition for resuspension
            #if TAU >= 0.05:
                # Resuspension probability (30% for marine particles, 70% for colloidal)
                if ParcelsRandom.uniform(0, 1) >= 0.3:
                    particle.status = 2  # Colloidal/Dissolved PBDEs
                else:
                    particle.status = 3  # Marine Particles         
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