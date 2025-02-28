#
### Ocean Parcels Kernels for PBDEs Simulations ###    
#
#
def PBDEs_states(particle, fieldset, time):
    #print('PBDEs_states kernel is running')
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
        if particle.time > (particle.time_release):
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
#### PBDEs states sinking velocities features ####
def PBDEs_forms(particle, fieldset, time):
    #print('PBDEs_forms kernel is running')
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
    # Note: U, V and e3t particle.depth are in meters, so I cannot use mbathy as an index, since this gives me Z levels, not depth!!
    # Since when particle.status == 4 the particle depth is set as the Total Depth. Therefore, if we do total depth - e3t we will get the depth of the second last grid cell...
    # which is the one that we want for the U and V calculations used in getting U_star!!
    if particle.status == 4:
    # If particle moves in depth, lon or lat, then update, if not, keep constant:
        if (particle.depth != particle.prev_depth or 
            particle.lat != particle.prev_lat or 
            particle.lon != particle.prev_lon):
            #              
            particle.e3t_val = fieldset.e3t[0, particle.depth, particle.lat, particle.lon]  # Update e3t if particle moves
            particle.log_e3t = math.log(particle.e3t_val * 0.5)  # Update log term if particle moves
            #    
            # Save previous location to keep it consistent
            particle.prev_depth = particle.depth
            particle.prev_lat = particle.lat
            particle.prev_lon = particle.lon
            #
        #  get the depth of the second last grid cell
        bat_particle = particle.depth - .5 *particle.e3t_val 
        #
        # Velocity in m/s
        u_vel = fieldset.U[time, bat_particle, particle.lat, particle.lon] * (particle.deg2met * particle.latT)
        particle.u_vel = u_vel
        v_vel = fieldset.V[time, bat_particle, particle.lat, particle.lon] * (particle.deg2met)
        particle.v_vel = v_vel
        # Horizontal velocity        
        # Algebraic arange for more efficiency
        TAU = ((u_vel ** 2) + (v_vel) ** 2) * (particle.k_constant) * ((particle.log_e3t - particle.log_z_star) ** (-2)) * 1024
        # Store Tau to see the numbers
        particle.tau_values = TAU
        # Apply resuspension 
        if TAU >= particle.threshold: # for colloids and marine particles
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