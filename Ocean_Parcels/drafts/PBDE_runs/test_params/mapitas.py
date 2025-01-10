from mpl_toolkits.basemap import Basemap
from matplotlib.transforms import Affine2D
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from cartopy import crs, feature
import zarr 
import xarray as xr
import numpy as np
#
#
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
#
clat = 49.195045
clon = -123.301956
#
######### NICE PLOTS ######### 
#
#
def plot_nice_map(data, figure, plot_axis, data_lon, data_lat, size, color_particle, label_particle, percentage):
    #
    rotation_angle = 10  # Rotation angle in degrees
    #
    #
    w_map = [-126.5, -121, 46.7, 51.3]
    # Define a rotated pole projection with a shifted pole to create the rotation effect
    m = Basemap(projection='lcc', resolution='h',
                lon_0=clon, lat_0=clat,
                llcrnrlon=w_map[0], urcrnrlon=w_map[1],
                llcrnrlat=w_map[2], urcrnrlat=w_map[3], ax=plot_axis)
    #
    lons = np.arange(*np.floor([w_map[0], w_map[1] + 1]))
    lats = np.arange(*np.floor([w_map[2], w_map[3] + 1]))
    #
    labels = [[0, 0, 1, 0], [1, 0, 0, 0]]
    #
    m.drawcoastlines(zorder=1)
    m.fillcontinents(color='Burlywood', zorder=0)
    m.drawmeridians(lons, color='k',labels = [0, 0, 0, 1], yoffset=None, zorder=2, fontsize = 12)
    m.drawparallels(lats, color='k',labels = labels[1], xoffset=None, zorder=2, fontsize = 12)
    m.drawrivers(zorder=2)
    #
    x, y = m(coords.nav_lon, coords.nav_lat)
    #
    blevels = list(np.arange(0,500,20))
    #
    C = plot_axis.contourf(x, y, mask.totaldepth[:,:], zorder=1,cmap = 'Greys',levels=blevels, extend='both')
    plot_axis.contourf(x, y, mask.totaldepth[:,:], [-0.01, 0.01], colors='lightgray', zorder=3)
    plot_axis.contour( x, y, mask.totaldepth[:,:], [0], colors='Black', zorder=4)
    #
    #C_lon, C_lat = m(clon, clat)
    plot_axis.plot(x[ :,  0], y[ :,  0], 'k-', zorder=6)
    plot_axis.plot(x[ :, -1], y[ :, -1], 'k-', zorder=6)
    plot_axis.plot(x[ 0,  :], y[ 0,  :], 'k-', zorder=6)
    plot_axis.plot(x[-1,  :], y[-1,  :], 'k-', zorder=6)
    #
    data_longitude, data_latitude = m(data_lon, data_lat)
    plot_axis.scatter(data_longitude, data_latitude, marker = '.', color = color_particle, s = size, label = label_particle)
    plot_axis.legend(loc='upper right', fontsize=16)
    #plot_axis.text(0.02, 0.01, str(np.round(percentage,2)) +' % ' + 'from Total', transform=plot_axis.transAxes, fontweight = 'bold', fontsize = 12)
    plt.suptitle('Spatial Distribution of PBDEs Phases between ' + str(data.time[0,0].values)[:10] + ' and ' + str(data.time[0,-1].values)[:10], fontsize = 18)
#
#
def plot_distribution(data_zarr, depth, status_number):
    #
    fig,axs=plt.subplots(1,2,figsize=(14,7))
    # Create a nice layout #
    rotation_angle = 10  # Rotation angle in degrees
    #
    #
    w_map = [-126.5, -121, 46.7, 51.3]
    # Define a rotated pole projection with a shifted pole to create the rotation effect
    m0 = Basemap(projection='lcc', resolution='h',
                lon_0=clon, lat_0=clat,
                llcrnrlon=w_map[0], urcrnrlon=w_map[1],
                llcrnrlat=w_map[2], urcrnrlat=w_map[3], ax=axs[0])
    #
    m1 = Basemap(projection='lcc', resolution='h',
                lon_0=clon, lat_0=clat,
                llcrnrlon=w_map[0], urcrnrlon=w_map[1],
                llcrnrlat=w_map[2], urcrnrlat=w_map[3], ax=axs[1])
    #    
    lons = np.arange(*np.floor([w_map[0], w_map[1] + 1]))
    lats = np.arange(*np.floor([w_map[2], w_map[3] + 1]))
    #
    labels = [[0, 0, 1, 0], [1, 0, 0, 0]]
    #
    m0.drawcoastlines(zorder=1)
    m0.fillcontinents(color='Burlywood', zorder=0)
    m0.drawmeridians(lons, color='k',labels = labels[0], yoffset=None, zorder=2)
    m0.drawparallels(lats, color='k',labels = labels[1], xoffset=None, zorder=2)
    m0.drawrivers(zorder=2)
    #
    m1.drawcoastlines(zorder=1)
    m1.fillcontinents(color='Burlywood', zorder=0)
    m1.drawmeridians(lons, color='k',labels = labels[0], yoffset=None, zorder=2)
    m1.drawparallels(lats, color='k',labels = labels[1], xoffset=None, zorder=2)
    m1.drawrivers(zorder=2)
    #    
    x, y = m0(coords.nav_lon, coords.nav_lat)
    axs[0].plot(x[ :,  0], y[ :,  0], 'k-', zorder=6)
    axs[0].plot(x[ :, -1], y[ :, -1], 'k-', zorder=6)
    axs[0].plot(x[ 0,  :], y[ 0,  :], 'k-', zorder=6)
    axs[0].plot(x[-1,  :], y[-1,  :], 'k-', zorder=6)
    #
    axs[1].plot(x[ :,  0], y[ :,  0], 'k-', zorder=6)
    axs[1].plot(x[ :, -1], y[ :, -1], 'k-', zorder=6)
    axs[1].plot(x[ 0,  :], y[ 0,  :], 'k-', zorder=6)
    axs[1].plot(x[-1,  :], y[-1,  :], 'k-', zorder=6)    
    #######################
    cmap = 'jet'
    #### Make map ####
    blevels = list(np.arange(0,500,20))
    #
    C = axs[0].contourf(x, y, mask.totaldepth[:,:], zorder=1,cmap = 'Greys',levels=blevels, extend='both')
    axs[0].contourf(x, y, mask.totaldepth[:,:], [-0.01, 0.01], colors='lightgray', zorder=3)
    axs[0].contour( x, y, mask.totaldepth[:,:], [0], colors='Black', zorder=4)
    #
    #
    data_lon0, data_lat0 = m0(data_zarr.lon.where(data_zarr.status == status_number), data_zarr.lat.where(data_zarr.status == status_number))
    #
    im = axs[0].scatter(data_lon0,data_lat0,zorder=3,c=depth.where(data_zarr.status == status_number), cmap = cmap ,s=2, vmin = 0, vmax = 500)
    #
    if status_number == 1:
        axs[0].text(0.02, 0.01, 'All Particles (S. Particle)', transform=axs[0].transAxes)
    elif status_number == 2:
        axs[0].text(0.02, 0.01, 'All Particles (Colloidal)', transform=axs[0].transAxes)        
    elif status_number == 3:
        axs[0].text(0.02, 0.01, 'All Particles (M. Particle)', transform=axs[0].transAxes)        
    else:
        axs[0].text(0.02, 0.01, 'All Particles (Bottom)', transform=axs[0].transAxes)        
    #
    #
    axs[1].contourf(x, y, mask.totaldepth[:,:], zorder=1,cmap = 'Greys',levels=blevels, extend='both')
    axs[1].contourf(x, y, mask.totaldepth[:,:], [-0.01, 0.01], colors='lightgray', zorder=3)
    axs[1].contour(x, y, mask.totaldepth[:,:], [0], colors='Black', zorder=4)
    #
    #
    data_lon1, data_lat1 = m1(data_zarr.lon[:,-1].where(data_zarr.status[:,-1] == status_number), data_zarr.lat[:,-1].where(data_zarr.status[:,-1] == status_number))
    #
    axs[1].scatter(data_lon1, data_lat1,zorder=3,c=depth[:,-1].where(data_zarr.status[:,-1] == status_number), cmap = cmap,s=2, vmin = 0, vmax = 500)
    #
    if status_number == 1:
        axs[1].text(0.02, 0.01, 'All Particles (S. Particle)', transform=axs[1].transAxes)        
    elif status_number == 2:
        axs[1].text(0.02, 0.01, 'Final Particles (Colloidal)', transform=axs[1].transAxes)        
    elif status_number == 3:
        axs[1].text(0.02, 0.01, 'Final Particles (M. Particle)', transform=axs[1].transAxes)
    else:
        axs[1].text(0.02, 0.01, 'Final Particles (Bottom)', transform=axs[1].transAxes)
    #
    cbar = fig.colorbar(im, ax=axs[1], location='right', shrink=0.8)
    cbar.set_label('Particles Depth [m]')
    #
    plt.suptitle('Spatial Distribution of PBDEs Phases between ' + str(data_zarr.time[0,0].values)[:10] + ' and ' + str(data_zarr.time[0,-1].values)[:10], fontsize = 15)
    plt.tight_layout()