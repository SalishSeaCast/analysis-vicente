#
### Ocean Parcels Kernels for PBDEs Simulations ###    
#
# Force kernels to not move particles when not released!! <---- Not to move or change state
#
def PBDEs_states(particle, fieldset, time):  

#*    particle.steps = 0
#*    if particle.status == particle.prestatus:
#*        particle.prestatus = particle.status + 0.10
#*    else:
#*        particle.prestatus = particle.status
    particle.H_vel_2 = 0
    particle.crit = 0
    particle.bat_particle = 0
    if (time > particle.release_time):
        if particle.status < 0:
            particle.status = - particle.status
        else:
#**            if particle.status > 100:
#**                particle.status -= 100
            random_value = ParcelsRandom.random()
        # Status updates
            if particle.status == 1 and random_value < fieldset.deso_s_probability:
                particle.status = 2  # Becomes Colloidal/Dissolved           
            elif particle.status == 11:
                if random_value < fieldset.sediment_burying_probability:
                    particle.status = 21
                elif random_value < fieldset.deso_s_probability:
                    particle.status = 12  # Becomes Colloidal/Dissolved in Sediments
            elif particle.status == 21 and random_value < fieldset.deso_s_probability:
                particle.status = 22  # Becomes Colloidal/Dissolved buried in sediments
                
            elif particle.status == 2 and random_value < fieldset.abso_probability:
                particle.status = 3  # Becomes attached to Marine Particle
            elif particle.status == 12:
                if random_value < fieldset.sediment_burying_probability:
                    particle.status = 22
                elif random_value < fieldset.abso_sed_probability:
                    particle.status = 13  # Becomes attached to Marine Particle in Sediments
            elif particle.status == 22 and random_value < fieldset.abso_sed_probability:
                particle.status = 23  # Becomes attached to Marine Particle buried in Sediments
                
            elif particle.status == 3 and random_value < fieldset.deso_m_probability:
                particle.status = 2  # Returns to Colloidal/Dissolved
            elif particle.status == 13:
                if random_value < fieldset.sediment_burying_probability:
                    particle.status = 23
                elif random_value < fieldset.deso_sed_probability:
                    particle.status = 12  # Returns to Colloidal/Dissolved in Sediments
            elif particle.status == 23 and random_value < fieldset.deso_sed_probability:
                particle.status = 22 # Returns to Colloidal/Dissolved buried Sediments
            
#### PBDEs states sinking velocities features ####
def Sinking(particle, fieldset, time):
    particle_ddepth = 0
    if particle.status == 1:
        particle_ddepth += fieldset.sinkvel_sewage * particle.dt
    #Sewage Particles sink fast        

    elif particle.status == 3:
        particle_ddepth += fieldset.sinkvel_marine * particle.dt

    # do settling here
    if particle.status == 1 or particle.status == 3:
        tda = fieldset.totaldepth[time, particle.depth, 
                        particle.lat, particle.lon]
        if particle_ddepth + particle.depth > tda:
            particle.depth  = tda # Get particles attached to the bottom when they reach it
            particle_ddepth = 0 # As I've put them on the bottom and that's where I want them.
            particle.status += 10 
    
#
#### ADVECTION ####
def Advection(particle, fieldset, time):
    #print('Advection kernel is running') 
    # Advection for all PBDEs in status 1, 2 and 3
    if particle.status == 1  or particle.status == 2 or particle.status == 3: 
        ssh = fieldset.sossheig[time, particle.depth, particle.lat, particle.lon] #SSH(t) sea surface height
        sshn = fieldset.sossheig[time+particle.dt, particle.depth, particle.lat, particle.lon] #SSH(t+dt) sea surface height in the next time step
        td = fieldset.totaldepth[time, particle.depth, particle.lat, particle.lon]#Total_depth 
        particle.fact = (1 + ssh / td)
        VVL = (sshn - ssh) * particle.depth / td # following Jose: 
        # https://github.com/SalishSeaCast/analysis-jose/blob/main/Graham/Source/OP_Kernels.py
        # calculate once and reuse
        dt_factor = 0.5 * particle.dt / particle.fact
        (u1, v1, w1) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
        lon1 = particle.lon + u1 * 0.5*particle.dt
        lat1 = particle.lat + v1 * 0.5*particle.dt
        dep1 = particle.depth + w1 * dt_factor
        #
        (u2, v2, w2) = fieldset.UVW[time + 0.5 * particle.dt, dep1, lat1, lon1]
        lon2 = particle.lon + u2 * 0.5*particle.dt
        lat2 = particle.lat + v2 * 0.5*particle.dt
        dep2 = particle.depth + w2 * dt_factor
        #
        (u3, v3, w3) = fieldset.UVW[time + 0.5 * particle.dt, dep2, lat2, lon2]
        lon3 = particle.lon + u3 * particle.dt
        lat3 = particle.lat + v3 * particle.dt
        dep3 = particle.depth + w3 * 2 * dt_factor
        #
        (u4, v4, w4) = fieldset.UVW[time + particle.dt, dep3, lat3, lon3]
        #
        wa = (w1 + 2*w2 + 2*w3 + w4) /6.
        particle_dlon = (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle_dlat = (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        particle_ddepth = particle_ddepth + wa/particle.fact * particle.dt + VVL
        
        if particle_ddepth + particle.depth < 0:
            particle_ddepth = - (2 * particle.depth + particle_ddepth)

#### TURBULENT MIX ####
def turb_mix(particle,fieldset,time):
    if particle.status == 1 or particle.status == 2 or particle.status == 3:
        """Vertical mixing"""
        #Vertical mixing
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
        # stick if mixed into bottom
        tdn = fieldset.totaldepth[time, particle.depth, 
                        particle.lat+particle_dlat, particle.lon+particle_dlon]
        if dzs + particle_ddepth + particle.depth > tdn:
            particle.depth  = tdn # Get particles attached to the bottom when they reach it
            particle_ddepth = 0 # As I've put them on the bottom and that's where I want them.
            particle.status += 10 
            particle.lat = particle.lat + particle_dlat
            particle_dlat = 0
            particle.lon = particle.lon + particle_dlon
            particle_dlon = 0
            #
        elif dzs + particle.depth + particle_ddepth < 0:
            particle_ddepth = -(dzs + 2*particle.depth + particle_ddepth) #reflection on surface
        #
        else:
            particle_ddepth += dzs #apply mixing
        
#
#### RESUSPENSION ####

def resuspension(particle, fieldset, time):
    # Only proceed if particle is at the seabed (status 11, 12 or 13, 21, 22, 23)
    if particle.status > 10:
        if particle.status > 20:
            vtau_constant = fieldset.tau_bury_constant
            vtau_constant_lower = fieldset.tau_bury_constant_lower
            vtau_constant_upper = fieldset.tau_bury_constant_upper
        else:
            vtau_constant = fieldset.tau_constant
            vtau_constant_lower = fieldset.tau_constant_lower
            vtau_constant_upper = fieldset.tau_constant_upper
        #
        
        tdr = fieldset.totaldepth[time, particle.depth, 
                        particle.lat, particle.lon]                 # even if new, already moved
        sshr = fieldset.sossheig[time, particle.depth, particle.lat, particle.lon] #SSH(t) sea surface height
        particle.fact = (1 + sshr / tdr)
        e3t_val_o2 = fieldset.e3t[time, tdr, particle.lat, particle.lon] * 0.5
        particle.e3t = e3t_val_o2
        bat_particle = max(tdr - e3t_val_o2, 0.5) 
        particle.bat_particle = bat_particle
        #
        # horizontal velocities in m/s  
        u_vel, v_vel = fieldset.UV[time, bat_particle, particle.lat, particle.lon]
        u_vel = u_vel * fieldset.u_deg2mps
        v_vel = v_vel * fieldset.v_deg2mps
        
        # squared horizontal velocity
        H_vel_2 = u_vel**2 + v_vel**2
        particle.H_vel_2 = H_vel_2
        suspend_rate = 2 * fieldset.sinkvel_marine # so they can actually resuspend
        if e3t_val_o2 < fieldset.lowere3t_o2:
            particle.crit = vtau_constant_lower
            if vtau_constant_lower <= H_vel_2:
                particle.status -= 10
                particle.stuck = 0
                particle.depth = tdr - min(e3t_val_o2, suspend_rate*particle.dt) / particle.fact 
            else:
                particle.stuck += 1
        elif e3t_val_o2 > fieldset.uppere3t_o2:
            particle.crit = vtau_constant_upper
            if vtau_constant_upper <= H_vel_2:
                particle.status -= 10
                particle.depth = tdr - min(e3t_val_o2, suspend_rate*particle.dt) / particle.fact
                particle.stuck = 0
            else:
                particle.stuck += 1
        else:
            log_e3t = math.log(e3t_val_o2 / particle.fact)
            particle.crit = vtau_constant * (log_e3t - fieldset.log_z_star) ** 2
            if vtau_constant * (log_e3t - fieldset.log_z_star) ** 2 <= H_vel_2:
                particle.status -= 10
                particle.depth = tdr - min(e3t_val_o2, suspend_rate*particle.dt) / particle.fact
                particle.stuck = 0
            else:
                particle.stuck += 1

    particle.tmask = fieldset.tmask[0, particle.depth, particle.lat, particle.lon]
    particle.umask = fieldset.umask[0, particle.depth, particle.lat, particle.lon]
    particle.vmask = fieldset.vmask[0, particle.depth, particle.lat, particle.lon]
    particle.fmask = fieldset.fmask[0, particle.depth, particle.lat, particle.lon]
    

    (uus, vvs, wws) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
    particle.uvalue = uus
    particle.vvalue = vvs
    particle.wvalue = wws

    particle.totaldepth = fieldset.totaldepth[time, particle.depth, particle.lat, particle.lon]
    if particle.status > 10:
        if (uus == 0) and (vvs == 0):
            if (particle.fmask == 0) or (particle.stuck > 86500/5*30):        
                dxp = fieldset.totaldepth[time, particle.depth, particle.lat+fieldset.dx_lat, particle.lon+fieldset.dx_lon]
                dxm = fieldset.totaldepth[time, particle.depth, particle.lat-fieldset.dx_lat, particle.lon-fieldset.dx_lon]
                dyp = fieldset.totaldepth[time, particle.depth, particle.lat+fieldset.dy_lat, particle.lon+fieldset.dy_lon]
                dym = fieldset.totaldepth[time, particle.depth, particle.lat-fieldset.dy_lat, particle.lon-fieldset.dy_lon]
    
                if dxp >= dyp:
                    if dxp >= dym:
                        if dxp >= dxm: # dxp is largest
                            doit = 0 # 'dxp_largest'
                        else: # dxm > dxp > dyp & dym
                            doit = 1 #'dxm_largest'
                    else: # dym > dxp > dyp
                        if dym >= dxm: #dym is largest
                            doit = 2 #'dym_largest'
                        else:
                            doit = 1 #'dxm_largest'
                else: # dyp > dxp
                    if dyp >= dym:
                        if dyp >= dxm:
                            doit = 3 #'dyp_largest'
                        else: #dxm > dyp > dxp & dym 
                            doit = 1 # 'dxm_largest'
                    else: #dym > dyp > dxp
                        if dym >= dxm:
                            doit = 2 #'dym_largest'
                        else:
                            doit = 1 #'dxm_largest'
                
                if doit == 0: #'dxp_largest':
                    particle.lat = particle.lat + dx_lat
                    particle.lon = particle.lon + dx_lon
                elif doit == 1:#'dxm_largest': 
                    particle.lat = particle.lat - dx_lat
                    particle.lon = particle.lon - dx_lon
                elif doit == 2:#'dym_largest':
                    particle.lat = particle.lat - dy_lat
                    particle.lon = particle.lon - dy_lon
                else:
                    particle.lat = particle.lat + dy_lat
                    particle.lon = particle.lon + dy_lon
                
                particle.totaldepth = fieldset.totaldepth[time, particle.depth, particle.lat, particle.lon]
                if particle.depth < particle.totaldepth:
                    particle.status -= 10
                    particle.stuck = 0
                
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

        