#
### Ocean Parcels Kernels for PBDEs Simulations ###    
#
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