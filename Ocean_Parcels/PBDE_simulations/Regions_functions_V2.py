#
import sys
import random
import xarray as xr
import numpy as np
import os
import yaml
import math
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import colors
import matplotlib.gridspec as gridspec
from cartopy import crs, feature
import zarr 
import glob
from matplotlib.colors import ListedColormap
# For polygons
from shapely.geometry import Polygon, Point
#
############ Grid X and Grid Y points for regions boundaries
# Northern Strait
# Subregion Left
Nx1_1, Ny1_1 = [115,185], [665,650]
Nx1_2, Ny1_2 = [110,220], [510,565]
#Subregion Right
Nx2_1, Ny2_1 = [228,258], [650,565]
Nx2_2, Ny2_2 = [220,258], [565,565]
Nx2_3, Ny2_3 = [185,228], [650,650]
###################################
# Central Strait
Cx1_1, Cy1_1 = [185,300], [505,505]
Cx1_2, Cy1_2 = [185,170], [505,540]
Cx1_3, Cy1_3 = [170,220], [540,565]
Cx1_4, Cy1_4 = [220,258], [565,565]
Cx1_5, Cy1_5 = [258,300], [565,505]
###################################
# Southern Strait
# Subregion North
Sx1_1, Sy1_1 = [202,292], [505,505]
Sx1_2, Sy1_2 = [300,320], [503,457]
Sx1_3, Sy1_3 = [210,238], [500,465]
Sx1_4, Sy1_4 = [238,245], [465,410]
Sx1_5, Sy1_5 = [320,302], [457,410]
Sx1_6, Sy1_6 = [245,302], [410,410]
# Subregion South
Sx2_1, Sy2_1 = [245,302], [410,410]
Sx2_2, Sy2_2 = [245,280], [410,350]
Sx2_3, Sy2_3 = [280,342], [350,315]
Sx2_4, Sy2_4 = [302,370], [410,390]
###################################
# JdF
Jx1_1, Jy1_1 = [145,145], [230,320]
Jx1_2, Jy1_2 = [145,80], [320,410]
Jx1_3, Jy1_3 = [145,60], [230,280]

# Puget Sound
Px1, Py1 = [203,312], [230,230]
# Haro Strait
Hx1_1, Hy1_1 = [200, 250],[305,300]
Hx1_2, Hy1_2 = [225, 242],[355,355]
Hx1_3, Hy1_3 = [244, 247],[355,320]
#
# Fraser River
Fx1_1, Fy1_1 = [320,302], [457,410]
Fx1_2, Fy1_2 = [302,370], [410,390]
Fx1_3, Fy1_3 = [320,370], [457,430]
#
# Howe Sound
HWx1_1, HWy1_1 = [300,320], [503,457]
HWx1_2, HWy1_2 = [300,350], [503,580]
HWx1_3, HWy1_3 = [350,395], [580,580]
HWx1_4, HWy1_4 = [395,395], [580,510]
HWx1_5, HWy1_5 = [395,345], [510,445]
HWx1_6, HWy1_6 = [345,320], [445,457]
#
###################################################
############ Paths #############
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
############## Longitude and Latitude boundary points ##############
# Northern Strait
# Subregion Left
N1_lon_1, N1_lat_1 = [mask.nav_lon[Ny1_1[0],Nx1_1[0]].values, mask.nav_lon[Ny1_1[1],Nx1_1[1]].values], [mask.nav_lat[Ny1_1[0],Nx1_1[0]].values, mask.nav_lat[Ny1_1[1],Nx1_1[1]].values]
N1_lon_2, N1_lat_2 = [mask.nav_lon[Ny1_2[0],Nx1_2[0]].values, mask.nav_lon[Ny1_2[1],Nx1_2[1]].values], [mask.nav_lat[Ny1_2[0],Nx1_2[0]].values, mask.nav_lat[Ny1_2[1],Nx1_2[1]].values]
#
# Subregion Right
N2_lon_1, N2_lat_1 = [mask.nav_lon[Ny2_1[0],Nx2_1[0]].values, mask.nav_lon[Ny2_1[1],Nx2_1[1]].values], [mask.nav_lat[Ny2_1[0],Nx2_1[0]].values, mask.nav_lat[Ny2_1[1],Nx2_1[1]].values]
N2_lon_2, N2_lat_2 = [mask.nav_lon[Ny2_2[0],Nx2_2[0]].values, mask.nav_lon[Ny2_2[1],Nx2_2[1]].values], [mask.nav_lat[Ny2_2[0],Nx2_2[0]].values, mask.nav_lat[Ny2_2[1],Nx2_2[1]].values]
N2_lon_3, N2_lat_3 = [mask.nav_lon[Ny2_3[0],Nx2_3[0]].values, mask.nav_lon[Ny2_3[1],Nx2_3[1]].values], [mask.nav_lat[Ny2_3[0],Nx2_3[0]].values, mask.nav_lat[Ny2_3[1],Nx2_3[1]].values]
#
# Central Strait
C1_lon_1, C1_lat_1 = [mask.nav_lon[Cy1_1[0],Cx1_1[0]].values, mask.nav_lon[Cy1_1[1],Cx1_1[1]].values], [mask.nav_lat[Cy1_1[0],Cx1_1[0]].values, mask.nav_lat[Cy1_1[1],Cx1_1[1]].values]
C1_lon_2, C1_lat_2 = [mask.nav_lon[Cy1_2[0],Cx1_2[0]].values, mask.nav_lon[Cy1_2[1],Cx1_2[1]].values], [mask.nav_lat[Cy1_2[0],Cx1_2[0]].values, mask.nav_lat[Cy1_2[1],Cx1_2[1]].values]
C1_lon_3, C1_lat_3 = [mask.nav_lon[Cy1_3[0],Cx1_3[0]].values, mask.nav_lon[Cy1_3[1],Cx1_3[1]].values], [mask.nav_lat[Cy1_3[0],Cx1_3[0]].values, mask.nav_lat[Cy1_3[1],Cx1_3[1]].values]
C1_lon_4, C1_lat_4 = [mask.nav_lon[Cy1_4[0],Cx1_4[0]].values, mask.nav_lon[Cy1_4[1],Cx1_4[1]].values], [mask.nav_lat[Cy1_4[0],Cx1_4[0]].values, mask.nav_lat[Cy1_4[1],Cx1_4[1]].values]
C1_lon_5, C1_lat_5 = [mask.nav_lon[Cy1_5[0],Cx1_5[0]].values, mask.nav_lon[Cy1_5[1],Cx1_5[1]].values], [mask.nav_lat[Cy1_5[0],Cx1_5[0]].values, mask.nav_lat[Cy1_5[1],Cx1_5[1]].values]
# Southern Strait
# Subregion North
S1_lon_1, S1_lat_1 = [mask.nav_lon[Sy1_1[0],Sx1_1[0]].values, mask.nav_lon[Sy1_1[1],Sx1_1[1]].values], [mask.nav_lat[Sy1_1[0],Sx1_1[0]].values, mask.nav_lat[Sy1_1[1],Sx1_1[1]].values]
S1_lon_2, S1_lat_2 = [mask.nav_lon[Sy1_2[0],Sx1_2[0]].values, mask.nav_lon[Sy1_2[1],Sx1_2[1]].values], [mask.nav_lat[Sy1_2[0],Sx1_2[0]].values, mask.nav_lat[Sy1_2[1],Sx1_2[1]].values]
S1_lon_3, S1_lat_3 = [mask.nav_lon[Sy1_3[0],Sx1_3[0]].values, mask.nav_lon[Sy1_3[1],Sx1_3[1]].values], [mask.nav_lat[Sy1_3[0],Sx1_3[0]].values, mask.nav_lat[Sy1_3[1],Sx1_3[1]].values]
S1_lon_4, S1_lat_4 = [mask.nav_lon[Sy1_4[0],Sx1_4[0]].values, mask.nav_lon[Sy1_4[1],Sx1_4[1]].values], [mask.nav_lat[Sy1_4[0],Sx1_4[0]].values, mask.nav_lat[Sy1_4[1],Sx1_4[1]].values]
S1_lon_5, S1_lat_5 = [mask.nav_lon[Sy1_5[0],Sx1_5[0]].values, mask.nav_lon[Sy1_5[1],Sx1_5[1]].values], [mask.nav_lat[Sy1_5[0],Sx1_5[0]].values, mask.nav_lat[Sy1_5[1],Sx1_5[1]].values]
S1_lon_6, S1_lat_6 = [mask.nav_lon[Sy1_6[0],Sx1_6[0]].values, mask.nav_lon[Sy1_6[1],Sx1_6[1]].values], [mask.nav_lat[Sy1_6[0],Sx1_6[0]].values, mask.nav_lat[Sy1_6[1],Sx1_6[1]].values]
# Subregion South
S2_lon_1, S2_lat_1 = [mask.nav_lon[Sy2_1[0],Sx2_1[0]].values, mask.nav_lon[Sy2_1[1],Sx2_1[1]].values], [mask.nav_lat[Sy2_1[0],Sx2_1[0]].values, mask.nav_lat[Sy2_1[1],Sx2_1[1]].values]
S2_lon_2, S2_lat_2 = [mask.nav_lon[Sy2_2[0],Sx2_2[0]].values, mask.nav_lon[Sy2_2[1],Sx2_2[1]].values], [mask.nav_lat[Sy2_2[0],Sx2_2[0]].values, mask.nav_lat[Sy2_2[1],Sx2_2[1]].values]
S2_lon_3, S2_lat_3 = [mask.nav_lon[Sy2_3[0],Sx2_3[0]].values, mask.nav_lon[Sy2_3[1],Sx2_3[1]].values], [mask.nav_lat[Sy2_3[0],Sx2_3[0]].values, mask.nav_lat[Sy2_3[1],Sx2_3[1]].values]
S2_lon_4, S2_lat_4 = [mask.nav_lon[Sy2_4[0],Sx2_4[0]].values, mask.nav_lon[Sy2_4[1],Sx2_4[1]].values], [mask.nav_lat[Sy2_4[0],Sx2_4[0]].values, mask.nav_lat[Sy2_4[1],Sx2_4[1]].values]
#
# Haro Strait
H1_lon_1, H1_lat_1 = [mask.nav_lon[Hy1_1[0],Hx1_1[0]].values, mask.nav_lon[Hy1_1[1],Hx1_1[1]].values], [mask.nav_lat[Hy1_1[0],Hx1_1[0]].values, mask.nav_lat[Hy1_1[1],Hx1_1[1]].values]
H1_lon_2, H1_lat_2 = [mask.nav_lon[Hy1_2[0],Hx1_2[0]].values, mask.nav_lon[Hy1_2[1],Hx1_2[1]].values], [mask.nav_lat[Hy1_2[0],Hx1_2[0]].values, mask.nav_lat[Hy1_2[1],Hx1_2[1]].values]
H1_lon_3, H1_lat_3 = [mask.nav_lon[Hy1_3[0],Hx1_3[0]].values, mask.nav_lon[Hy1_3[1],Hx1_3[1]].values], [mask.nav_lat[Hy1_3[0],Hx1_3[0]].values, mask.nav_lat[Hy1_3[1],Hx1_3[1]].values]#
#
# Juan de Fuca Strait
J1_lon_1, J1_lat_1 = [mask.nav_lon[Jy1_1[0],Jx1_1[0]].values, mask.nav_lon[Jy1_1[1],Jx1_1[1]].values], [mask.nav_lat[Jy1_1[0],Jx1_1[0]].values, mask.nav_lat[Jy1_1[1],Jx1_1[1]].values]
J1_lon_2, J1_lat_2 = [mask.nav_lon[Jy1_2[0],Jx1_2[0]].values, mask.nav_lon[Jy1_2[1],Jx1_2[1]].values], [mask.nav_lat[Jy1_2[0],Jx1_2[0]].values, mask.nav_lat[Jy1_2[1],Jx1_2[1]].values]
J1_lon_3, J1_lat_3 = [mask.nav_lon[Jy1_3[0],Jx1_3[0]].values, mask.nav_lon[Jy1_3[1],Jx1_3[1]].values], [mask.nav_lat[Jy1_3[0],Jx1_3[0]].values, mask.nav_lat[Jy1_3[1],Jx1_3[1]].values]#
#
# Fraser River
F1_lon_1, F1_lat_1 = [mask.nav_lon[Fy1_1[0],Fx1_1[0]].values, mask.nav_lon[Fy1_1[1],Fx1_1[1]].values], [mask.nav_lat[Fy1_1[0],Fx1_1[0]].values, mask.nav_lat[Fy1_1[1],Fx1_1[1]].values]
F1_lon_2, F1_lat_2 = [mask.nav_lon[Fy1_2[0],Fx1_2[0]].values, mask.nav_lon[Fy1_2[1],Fx1_2[1]].values], [mask.nav_lat[Fy1_2[0],Fx1_2[0]].values, mask.nav_lat[Fy1_2[1],Fx1_2[1]].values]
F1_lon_3, F1_lat_3 = [mask.nav_lon[Fy1_3[0],Fx1_3[0]].values, mask.nav_lon[Fy1_3[1],Fx1_3[1]].values], [mask.nav_lat[Fy1_3[0],Fx1_3[0]].values, mask.nav_lat[Fy1_3[1],Fx1_3[1]].values]#
#
# Howe Sound
HW1_lon_1, HW1_lat_1 = [mask.nav_lon[HWy1_1[0],HWx1_1[0]].values, mask.nav_lon[HWy1_1[1],HWx1_1[1]].values], [mask.nav_lat[HWy1_1[0],HWx1_1[0]].values, mask.nav_lat[HWy1_1[1],HWx1_1[1]].values]
HW1_lon_2, HW1_lat_2 = [mask.nav_lon[HWy1_2[0],HWx1_2[0]].values, mask.nav_lon[HWy1_2[1],HWx1_2[1]].values], [mask.nav_lat[HWy1_2[0],HWx1_2[0]].values, mask.nav_lat[HWy1_2[1],HWx1_2[1]].values]
HW1_lon_3, HW1_lat_3 = [mask.nav_lon[HWy1_3[0],HWx1_3[0]].values, mask.nav_lon[HWy1_3[1],HWx1_3[1]].values], [mask.nav_lat[HWy1_3[0],HWx1_3[0]].values, mask.nav_lat[HWy1_3[1],HWx1_3[1]].values]
HW1_lon_4, HW1_lat_4 = [mask.nav_lon[HWy1_4[0],HWx1_4[0]].values, mask.nav_lon[HWy1_4[1],HWx1_4[1]].values], [mask.nav_lat[HWy1_4[0],HWx1_4[0]].values, mask.nav_lat[HWy1_4[1],HWx1_4[1]].values]
HW1_lon_5, HW1_lat_5 = [mask.nav_lon[HWy1_5[0],HWx1_5[0]].values, mask.nav_lon[HWy1_5[1],HWx1_5[1]].values], [mask.nav_lat[HWy1_5[0],HWx1_5[0]].values, mask.nav_lat[HWy1_5[1],HWx1_5[1]].values]
HW1_lon_6, HW1_lat_6 = [mask.nav_lon[HWy1_6[0],HWx1_6[0]].values, mask.nav_lon[HWy1_6[1],HWx1_6[1]].values], [mask.nav_lat[HWy1_6[0],HWx1_6[0]].values, mask.nav_lat[HWy1_6[1],HWx1_6[1]].values]

############# POLYGONS PER REGION ##########
#
polygon_lon_lat_N1 = [
    (N1_lon_1[0], N1_lat_1[0]),
    (N1_lon_1[1], N1_lat_1[1]),
    (N1_lon_2[1], N1_lat_2[1]),
    (N1_lon_2[0], N1_lat_2[0]),
    (N1_lon_1[0], N1_lat_1[0])
]
polygon_coors_N1 = Polygon(polygon_lon_lat_N1)
#Subregion Right
polygon_lon_lat_N2 = [
    (N2_lon_3[0], N2_lat_3[0]),
    (N2_lon_1[0], N2_lat_1[0]),
    (N2_lon_2[1], N2_lat_2[1]),
    (N2_lon_2[0], N2_lat_2[0])
]
polygon_coors_N2 = Polygon(polygon_lon_lat_N2)
#**Central Strait Polygon**
polygon_lon_lat_C1 = [
    (C1_lon_1[1], C1_lat_1[1]),
    (C1_lon_1[0], C1_lat_1[0]),
    (C1_lon_2[0], C1_lat_2[0]),
    (C1_lon_2[1], C1_lat_2[1]),
    (C1_lon_3[1], C1_lat_3[1]),
    (C1_lon_4[1], C1_lat_4[1])
]
polygon_coors_C1 = Polygon(polygon_lon_lat_C1)
#**Southern Strait Polygon**
#Subregion North
polygon_lon_lat_S1 = [
    (S1_lon_1[0], S1_lat_1[0]),
    (S1_lon_1[1], S1_lat_1[1]),
    (S1_lon_2[1], S1_lat_2[1]),
    (S1_lon_5[1], S1_lat_5[1]),
    (S1_lon_6[0], S1_lat_6[0]),
    (S1_lon_4[0], S1_lat_4[0])
]
polygon_coors_S1 = Polygon(polygon_lon_lat_S1)
#Subregion South
polygon_lon_lat_S2 = [
    (S2_lon_1[1], S2_lat_1[1]),
    (S2_lon_1[0], S2_lat_1[0]),
    (S2_lon_2[1], S2_lat_2[1]),
    (S2_lon_3[1], S2_lat_3[1]),
    (S2_lon_4[1], S2_lat_4[1])
]
polygon_coors_S2 = Polygon(polygon_lon_lat_S2)
#**Haro Strait Polygon**
polygon_lon_lat_H1 = [
    (H1_lon_1[1], H1_lat_1[1]),
    (H1_lon_1[0], H1_lat_1[0]),
    (H1_lon_2[0], H1_lat_2[0]),
    (H1_lon_3[0], H1_lat_3[0])

]
polygon_coors_H1 = Polygon(polygon_lon_lat_H1)
#**Juan de Fuca Strait Polygon**
polygon_lon_lat_J1 = [
    (J1_lon_1[0], J1_lat_1[0]),
    (J1_lon_1[1], J1_lat_1[1]),
    (J1_lon_2[1], J1_lat_2[1]),
    (J1_lon_3[1], J1_lat_3[1])

]
polygon_coors_J1 = Polygon(polygon_lon_lat_J1)
#**Fraser River**
polygon_lon_lat_F1 = [
    (F1_lon_1[0], F1_lat_1[0]),
    (F1_lon_1[1], F1_lat_1[1]),
    (F1_lon_2[1], F1_lat_2[1]),
    (F1_lon_3[1], F1_lat_3[1])
]
polygon_coors_F1 = Polygon(polygon_lon_lat_F1)
#**Howe Sound**
polygon_lon_lat_HW1 = [
    (HW1_lon_1[1], HW1_lat_1[1]),
    (HW1_lon_1[0], HW1_lat_1[0]),
    (HW1_lon_2[1], HW1_lat_2[1]),
    (HW1_lon_3[1], HW1_lat_3[1]),
    (HW1_lon_4[1], HW1_lat_4[1]),
    (HW1_lon_5[1], HW1_lat_5[1]),
    
]
polygon_coors_HW1 = Polygon(polygon_lon_lat_HW1)
#
def inside_polygon_lon_lat(polygon, mask_path = path['mask']):
    #
    mask = xr.open_dataset(mask_path)
    nav_lon = mask['nav_lon'].values
    nav_lat = mask['nav_lat'].values
    tmask = mask['tmask'][0, 0].values
    # 
    nav_lon_flat = nav_lon.ravel()
    nav_lat_flat = nav_lat.ravel()
    tmask_flat = tmask.ravel()
    #
    points = np.array([Point(lon, lat) for lon, lat in zip(nav_lon_flat, nav_lat_flat)])
    #
    inside_mask = np.array([polygon.contains(point) for point in points]) & (tmask_flat == 1)
    #
    longitudes = nav_lon_flat[inside_mask]
    latitudes = nav_lat_flat[inside_mask]
    #
    return longitudes, latitudes 
#
#
lon_NSoG_N1, lat_NSoG_N1 = inside_polygon_lon_lat(polygon_coors_N1)
lon_NSoG_N2, lat_NSoG_N2 = inside_polygon_lon_lat(polygon_coors_N2)
lon_CSoG_C1, lat_CSoG_C1 = inside_polygon_lon_lat(polygon_coors_C1)
lon_SSoG_S1, lat_SSoG_S1 = inside_polygon_lon_lat(polygon_coors_S1)
lon_SSoG_S2, lat_SSoG_S2 = inside_polygon_lon_lat(polygon_coors_S2)
lon_Haro_H1, lat_Haro_H1 = inside_polygon_lon_lat(polygon_coors_H1)
lon_Juan_J1, lat_Juan_J1 = inside_polygon_lon_lat(polygon_coors_J1)
lon_Fraser_F1, lat_Fraser_F1 = inside_polygon_lon_lat(polygon_coors_F1)
lon_Howe_HW1, lat_Howe_HW1 = inside_polygon_lon_lat(polygon_coors_HW1)
#
######### INPUT FILE FUNCTION ##########
def input_file(filename):
    data = xr.open_zarr(filename)
    return data
######### points inside each polygon ########
def points_inside(polygon, data, t):
    #
    lon_values = data.lon[:, t].values
    lat_values = data.lat[:, t].values
    depth_values =  data.z[:,t].values / data.fact[:,t].values
    status_values = data.status[:, t].values
    #
    points = np.array([Point(lon, lat) for lon, lat in zip(lon_values, lat_values)])
    #    
    inside_mask = np.array([polygon.contains(point) for point in points])
    #
    longitudes = lon_values[inside_mask]
    latitudes = lat_values[inside_mask]
    status_inside = status_values[inside_mask].astype(int)
    depth_inside = depth_values[inside_mask]
    #
    amount = np.sum(inside_mask)
    #
    return depth_inside, longitudes, latitudes, status_inside, amount
#
def polygon_definition(filename):
    # Function updated to get the info of the last day for each month
    data = input_file(filename)

    time_values = pd.to_datetime(data.time[0, :])
    
    df_time = pd.DataFrame({'time': time_values})
    df_time['month'] = df_time['time'].dt.to_period('M')
    last_of_month_indices = df_time.groupby('month').tail(1).index

    outputs = ['depth', 'lon', 'lat', 'status', 'n_particles']
    n_steps = len(last_of_month_indices)

    polygons_dict = {
        'N1': (polygon_coors_N1, pd.DataFrame(columns=outputs, index=np.arange(0, n_steps), dtype=object)),
        'N2': (polygon_coors_N2, pd.DataFrame(columns=outputs, index=np.arange(0, n_steps), dtype=object)),
        'C1': (polygon_coors_C1, pd.DataFrame(columns=outputs, index=np.arange(0, n_steps), dtype=object)),
        'S1': (polygon_coors_S1, pd.DataFrame(columns=outputs, index=np.arange(0, n_steps), dtype=object)),
        'S2': (polygon_coors_S2, pd.DataFrame(columns=outputs, index=np.arange(0, n_steps), dtype=object)),
        'H1': (polygon_coors_H1, pd.DataFrame(columns=outputs, index=np.arange(0, n_steps), dtype=object)),
        'J1': (polygon_coors_J1, pd.DataFrame(columns=outputs, index=np.arange(0, n_steps), dtype=object)),
        'F1': (polygon_coors_F1, pd.DataFrame(columns=outputs, index=np.arange(0, n_steps), dtype=object)),
        'HW1': (polygon_coors_HW1, pd.DataFrame(columns=outputs, index=np.arange(0, n_steps), dtype=object))
    }

    # Loop through only the end-of-month steps
    for step_index, i in enumerate(last_of_month_indices):
        for key, (polygon, data_sets) in polygons_dict.items():
            depth_region, lon_region, lat_region, status_region, n_part = points_inside(polygon, data, t=i)
            data_sets.iloc[step_index] = [depth_region, lon_region, lat_region, status_region, n_part]

    return polygons_dict

#
def vertical_mean_total_profiles(polygon_section, v_resolution):
    PROF = np.linspace(0, int(mask['totaldepth'].max().values), v_resolution)

    # DataFrames to store average and total profiles
    DATA_depth_mean = pd.DataFrame(columns=['Avg. Depth', 'Avg. Particles'])
    DATA_depth_total = pd.DataFrame(columns=['Avg. Depth', 'Total Particles'])

    for k in range(len(PROF)-1):
        bin_depths = []
        particles_per_row = []
        total_particles_in_bin = 0

        for jj in range(len(polygon_section[1])):
            row = polygon_section[1].iloc[jj]

            if len(row.depth) == 0:
                continue

            depths = np.array(row.depth)
            statuses = np.array(row.status)

            valid_mask = statuses > 0
            depths = depths[valid_mask]

            depths_in_bin = depths[(depths >= PROF[k]) & (depths < PROF[k+1])]

            if len(depths_in_bin) > 0:
                bin_depths.extend(depths_in_bin)
                particles_per_row.append(len(depths_in_bin))
                total_particles_in_bin += len(depths_in_bin)

        avg_depth = np.nanmean(bin_depths) if bin_depths else np.nan
        avg_particles = np.nanmean(particles_per_row) if particles_per_row else 0

        DATA_depth_mean.at[k, 'Avg. Depth'] = avg_depth
        DATA_depth_mean.at[k, 'Avg. Particles'] = avg_particles
        DATA_depth_total.at[k, 'Avg. Depth'] = avg_depth
        DATA_depth_total.at[k, 'Total Particles'] = total_particles_in_bin

    # Handle last bin
    k = len(PROF)-1
    bin_depths = []
    particles_per_row = []
    total_particles_in_bin = 0

    for jj in range(len(polygon_section[1])):
        row = polygon_section[1].iloc[jj]

        if len(row.depth) == 0:
            continue

        depths = np.array(row.depth)
        statuses = np.array(row.status)

        valid_mask = statuses > 0
        depths = depths[valid_mask]

        depths_in_bin = depths[(depths >= PROF[-2]) & (depths < PROF[-1])]

        if len(depths_in_bin) > 0:
            bin_depths.extend(depths_in_bin)
            particles_per_row.append(len(depths_in_bin))
            total_particles_in_bin += len(depths_in_bin)

    avg_depth = np.nanmean(bin_depths) if bin_depths else np.nan
    avg_particles = np.nanmean(particles_per_row) if particles_per_row else 0

    DATA_depth_mean.at[k, 'Avg. Depth'] = avg_depth
    DATA_depth_mean.at[k, 'Avg. Particles'] = avg_particles
    DATA_depth_total.at[k, 'Avg. Depth'] = avg_depth
    DATA_depth_total.at[k, 'Total Particles'] = total_particles_in_bin

    return DATA_depth_mean, DATA_depth_total

clat = [49.195045]
clon = [-123.301956]
#
nav_lon = mask['nav_lon'].values
nav_lat = mask['nav_lat'].values
tmask = mask['tmask'][0, 0].values
#
#
def vertical_status_profiles(polygon_section, v_resolution):

    PROF = np.linspace(0, int(mask['totaldepth'].max().values), v_resolution)

    status_list = [1, 2, 3, 11, 12, 13, 21, 22, 23]

    # Initialize the result DataFrame
    DATA_depth_mean = pd.DataFrame(
        columns=['Avg. Depth'] + [f'Particles Status {s}' for s in status_list],
        index=range(len(PROF) - 1)
    )
    DATA_depth_mean[:] = np.nan

    # Loop through depth bins
    for k in range(len(PROF) - 1):
        bin_depths = {s: [] for s in status_list}
        particles_per_status = {s: 0 for s in status_list}

        for i in range(len(polygon_section[1])):
            particle = polygon_section[1].iloc[i]
            statuses = particle.status
            depths = particle.depth

            if isinstance(statuses, (np.ndarray, list)) and isinstance(depths, (np.ndarray, list)) and len(statuses) == len(depths):
                for status, depth in zip(statuses, depths):
                    if status in status_list and PROF[k] <= depth < PROF[k + 1]:
                        bin_depths[status].append(depth)
                        particles_per_status[status] += 1

        #
        valid_depths = [np.nanmean(dlist) for dlist in bin_depths.values() if dlist]
        avg_depth = np.nanmean(valid_depths) if valid_depths else np.nan
        DATA_depth_mean.at[k, 'Avg. Depth'] = avg_depth

        for status in status_list:
            DATA_depth_mean.at[k, f'Particles Status {status}'] = particles_per_status[status]

    return DATA_depth_mean
     
#
def plot_vertical_total_profiles(array_vertical_profiles, data):
    #
    cmap = ListedColormap(['grey', 'white'])
    #
    fig = plt.figure(figsize=(15, 18))
    gs = fig.add_gridspec(5, 3, width_ratios=[1.0, 1.0, 0.8])
    #
    ax_plots = [fig.add_subplot(gs[i // 2, i % 2]) for i in range(9)]
    #
    ax_map = fig.add_subplot(gs[1:3, 2]) 
    #
    ax_plots[0].plot(array_vertical_profiles[0]['Total Particles'], array_vertical_profiles[0]['Avg. Depth'], '.-r', label='NSoG N1')
    ax_plots[1].plot(array_vertical_profiles[1]['Total Particles'], array_vertical_profiles[1]['Avg. Depth'], '.-b', label='NSoG N2')
    ax_plots[2].plot(array_vertical_profiles[2]['Total Particles'], array_vertical_profiles[2]['Avg. Depth'], '.-c', label='CSoG C1')
    ax_plots[3].plot(array_vertical_profiles[3]['Total Particles'], array_vertical_profiles[3]['Avg. Depth'], '.-g', label='SSoG S1')
    ax_plots[4].plot(array_vertical_profiles[4]['Total Particles'], array_vertical_profiles[4]['Avg. Depth'], 'tab:orange',linestyle='-', marker='.', label='SSoG S2')
    ax_plots[5].plot(array_vertical_profiles[5]['Total Particles'], array_vertical_profiles[5]['Avg. Depth'], '.-y', label='Haro H1')
    ax_plots[6].plot(array_vertical_profiles[6]['Total Particles'], array_vertical_profiles[6]['Avg. Depth'], '.-', color = 'm', label='JdF J1')
    ax_plots[7].plot(array_vertical_profiles[7]['Total Particles'], array_vertical_profiles[7]['Avg. Depth'], '.-', color = 'tab:brown', label='Fraser F1')
    ax_plots[8].plot(array_vertical_profiles[8]['Total Particles'], array_vertical_profiles[8]['Avg. Depth'], '.-', color = 'tab:gray', label='Howe Sound HW1')
    
    #
    for ax in ax_plots:
        ax.set_ylim(0, 430)
        ax.invert_yaxis()
        ax.legend(fontsize=10)
        ax.set_xlabel("Total N of Particles")
        ax.set_ylabel("Avg. Depth bins [m]")
        ax.grid(linestyle = '--')
    #
    ax_map.pcolormesh(nav_lon, nav_lat, tmask, cmap=cmap, shading='auto')
    ax_map.scatter(lon_NSoG_N1, lat_NSoG_N1, s=10, color='r', alpha=0.5, label='NSoG N1')
    ax_map.scatter(lon_NSoG_N2, lat_NSoG_N2, s=10, color='b', alpha=0.5, label='NSoG N2')
    ax_map.scatter(lon_CSoG_C1, lat_CSoG_C1, s=10, color='c', alpha=0.5, label='CSoG')
    ax_map.scatter(lon_SSoG_S1, lat_SSoG_S1, s=10, color='g', alpha=0.5, label='SSoG S1')
    ax_map.scatter(lon_SSoG_S2, lat_SSoG_S2, s=10, color='tab:orange', alpha=0.5, label='SSoG S2')
    ax_map.scatter(lon_Haro_H1, lat_Haro_H1, s=10, color='y', alpha=0.5, label='Haro H1')
    ax_map.scatter(lon_Juan_J1, lat_Juan_J1, s=10, color='m', alpha=0.5, label='JdF J1')
    ax_map.scatter(lon_Fraser_F1, lat_Fraser_F1, s=10, color='tab:brown', alpha=0.5, label='Fraser F1')
    ax_map.scatter(lon_Howe_HW1, lat_Howe_HW1, s=10, color='tab:gray', alpha=0.5, label='Howe Sound HW1')
    ax_map.scatter(data.lon, data.lat, s=.2, color='k', alpha = .2)
    ax_map.scatter(clon[0], clat[0], s=10, color='m', marker='s')
    #
    ax_map.set_xlim(-125.3, -122.2)
    ax_map.set_ylim(47.5, 50.5)
    ax_map.legend(fontsize=12)
    ax_map.tick_params(
        which='both', bottom=False, top=False, left=False, right=False,
        labelbottom=False, labelleft=False,
    )
    #
    fig.tight_layout()
    plt.show()      
#

def plot_vertical_status_profiles(status_label, status_profiles, states, colors, regions, data):
    fig = plt.figure(figsize=(18, 14))
    gs = fig.add_gridspec(3, 4, width_ratios=[1.0, 1.0, 1.0, 1.2])

    ax_plots = [fig.add_subplot(gs[i, j]) for i in range(3) for j in range(3)]

    ax_map = fig.add_subplot(gs[:, 3])

    for i, status_name in enumerate(status_label):
        ax = ax_plots[i]  
        for profile in range(len(status_profiles)):
            df = status_profiles[profile]

            if status_name not in df.columns:
                # if data is missing
                print(f"Warning: '{status_name}' not in profile {regions[profile]}")
                continue

            # Use twinx for S1
            if profile == 3:
                twin_axis = ax.twiny()
                twin_axis.plot(df[status_name], df['Avg. Depth'], colors[profile],
                               label=regions[profile], marker='.')
                twin_axis.tick_params(axis='x', colors='green')
                twin_axis.legend(loc='upper right')
            else:
                ax.plot(df[status_name], df['Avg. Depth'], colors[profile],
                        label=regions[profile], marker='.')

        ax.set_title(states[i])
        ax.set_ylim(0, 430)
        ax.invert_yaxis()
        ax.set_xlabel("Particle Count")
        ax.set_ylabel("Avg. Depth [m]")
        ax.grid(linestyle='--')
        ax.legend(loc='lower right', fontsize=9)

    # Plot map
    ax_map.pcolormesh(nav_lon, nav_lat, tmask, cmap=ListedColormap(['grey', 'white']), shading='auto')
    ax_map.scatter(lon_NSoG_N1, lat_NSoG_N1, s=10, color='r', alpha=0.5, label='NSoG (N1)')
    ax_map.scatter(lon_NSoG_N2, lat_NSoG_N2, s=10, color='b', alpha=0.5, label='NSoG (N2)')
    ax_map.scatter(lon_CSoG_C1, lat_CSoG_C1, s=10, color='c', alpha=0.5, label='CSoG (C1)')
    ax_map.scatter(lon_SSoG_S1, lat_SSoG_S1, s=10, color='g', alpha=0.5, label='SSoG (S1)')
    ax_map.scatter(lon_SSoG_S2, lat_SSoG_S2, s=10, color='tab:orange', alpha=0.5, label='SSoG (S2)')
    ax_map.scatter(lon_Haro_H1, lat_Haro_H1, s=10, color='y', alpha=0.5, label='Haro (H1)')
    ax_map.scatter(lon_Juan_J1, lat_Juan_J1, s=10, color='m', alpha=0.5, label='JdF (J1)')
    ax_map.scatter(lon_Fraser_F1, lat_Fraser_F1, s=10, color='tab:brown', alpha=0.5, label='Fraser (F1)')
    ax_map.scatter(lon_Howe_HW1, lat_Howe_HW1, s=10, color='tab:gray', alpha=0.5, label='Howe Sound (HW1)')    
    ax_map.scatter(data.lon, data.lat, s=.2, color='k', alpha = .2)
    ax_map.scatter(clon[0], clat[0], s=10, color='m', marker='s')

    ax_map.set_xlim(-125.3, -122.2)
    ax_map.set_ylim(47.5, 50.5)
    ax_map.legend(fontsize=11)
    ax_map.tick_params(which='both', bottom=False, top=False, left=False, right=False,
                       labelbottom=False, labelleft=False)

    fig.tight_layout()
    plt.show()
    #
    #
#

def plot_vertical_total_state_profiles(array_total_profiles, array_status_profiles,
                                       water_regions_lon, water_regions_lat,
                                       sedimented_regions_lon, sedimented_regions_lat,
                                       water_regions_depth, sedimented_Regions_depth):
    plt.rcParams.update({'font.size': 12})
    cmap = ListedColormap(['grey', 'white'])

    site_names = ['N1', 'N2', 'C1', 'S1', 'S2',
                  'H1', 'J1', 'F1', 'HW1']

    site_coords = [
        (lon_NSoG_N1, lat_NSoG_N1),
        (lon_NSoG_N2, lat_NSoG_N2),
        (lon_CSoG_C1, lat_CSoG_C1),
        (lon_SSoG_S1, lat_SSoG_S1),
        (lon_SSoG_S2, lat_SSoG_S2),
        (lon_Haro_H1, lat_Haro_H1),
        (lon_Juan_J1, lat_Juan_J1),
        (lon_Fraser_F1, lat_Fraser_F1),
        (lon_Howe_HW1, lat_Howe_HW1),
    ]

    fig, axes = plt.subplots(nrows=9, ncols=2, figsize=(12, 45), gridspec_kw={'width_ratios': [2, 1]})

    for i in range(9):
        ax_plot = axes[i, 0]
        ax_map = axes[i, 1]

        df_total = array_total_profiles[i]
        df_status = array_status_profiles[i]

        # Convert depth and particle count data
        y = df_total['Avg. Depth'].astype(float).to_numpy()
        y_status = df_status['Avg. Depth'].astype(float).to_numpy()
        total = df_total['Total Particles'].astype(float).to_numpy()
        water_col = (
            df_status['Particles Status 1'].astype(float) +
            df_status['Particles Status 2'].astype(float) +
            df_status['Particles Status 3'].astype(float)
        ).to_numpy()
        sedimented = (
            df_status['Particles Status 11'].astype(float) +
            df_status['Particles Status 12'].astype(float) +
            df_status['Particles Status 13'].astype(float)
        ).to_numpy()

        # Left: vertical profile
        ax_plot.fill_betweenx(y_status, water_col, color='blue', alpha=0.6, label='Water Column')
        ax_plot.fill_betweenx(y_status, sedimented, color='red', alpha=0.6, label='Sedimented')
        ax_plot.plot(total, y, color='k', linewidth=2, label='Total')
        ax_plot.set_ylim(0, 430)
        ax_plot.invert_yaxis()
        ax_plot.set_title(rf"$\bf{{{site_names[i]}}}$")
        ax_plot.set_xlabel("Particle Count")
        ax_plot.set_ylabel("Avg. Depth [m]")
        ax_plot.grid(True, linestyle='--')
        ax_plot.legend(fontsize=8, loc='lower right')

        # Right: map
        ax_map.pcolormesh(nav_lon, nav_lat, tmask, cmap=cmap, shading='auto')
        ax_map.scatter(*site_coords[i], s=15, alpha=0.03, color='k', label=site_names[i])
        ax_map.scatter(water_regions_lon[i], water_regions_lat[i], s=0.5, color='blue', alpha=0.5, label='Water')
        ax_map.scatter(sedimented_regions_lon[i], sedimented_regions_lat[i], s=0.5, color='red', alpha=0.5, label='Sedimented')
        ax_map.set_xlim(site_coords[i][0].mean() - 0.7, site_coords[i][0].mean() + 0.7)
        ax_map.set_ylim(site_coords[i][1].mean() - 0.7, site_coords[i][1].mean() + 0.7)
        ax_map.set_aspect('equal', adjustable='box')
        ax_map.set_title(rf"$\bf{{{site_names[i]}}}$")

        ax_map.tick_params(
            which='both', bottom=False, top=False, left=False, right=False,
            labelbottom=False, labelleft=False,
        )

        # Add depth info text box
        #mean_water_depth = np.nanmean(water_regions_depth[i])
        #mean_sediment_depth = np.nanmean(sedimented_Regions_depth[i])
        #mean_total_depth = np.nanmean(np.concatenate([water_regions_depth[i], sedimented_Regions_depth[i]]))

        #depth_text = (
        #    f"Mean Depth Total: {mean_total_depth:.1f} m\n"
        #    f"Mean Depth Water: {mean_water_depth:.1f} m\n"
        #    f"Mean Depth Sediment: {mean_sediment_depth:.1f} m"
        #)

        #ax_map.text(0.5, 0.02, depth_text, fontsize=9, color='k',
        #            ha='center', va='bottom', transform=ax_map.transAxes,
        #            bbox=dict(boxstyle="round,pad=0.3", facecolor='white', edgecolor='gray', alpha=0.8))

    fig.tight_layout()
    plt.show()


#
# Plot for Total particles, state and status ### :D
#
def plot_vertical_state_status(array_total_profiles, array_status_profiles):
    plt.rcParams.update({'font.size': 16})

    fig, axs = plt.subplots(nrows=9, ncols=3, figsize=(20, 40), sharey=True)
    axs = axs.reshape(9, 3)

    site_names = ['N1', 'N2', 'C1', 'S1', 'S2', 'H1', 'J1', 'F1', 'HW1']

    for i, site in enumerate(site_names):
        df_total = array_total_profiles[i]
        df_status = array_status_profiles[i]

        y = df_total['Avg. Depth'].astype(float).to_numpy()
        y_status = df_status['Avg. Depth'].astype(float).to_numpy()
        total = df_total['Total Particles'].astype(float).to_numpy()

        # Water column statuses
        status1 = df_status['Particles Status 1'].astype(float).to_numpy()
        status2 = df_status['Particles Status 2'].astype(float).to_numpy()
        status3 = df_status['Particles Status 3'].astype(float).to_numpy()
        water_col = status1 + status2 + status3

        # Sedimented statuses
        status11 = df_status['Particles Status 11'].astype(float).to_numpy()
        status12 = df_status['Particles Status 12'].astype(float).to_numpy()
        status13 = df_status['Particles Status 13'].astype(float).to_numpy()
        sedimented = status11 + status12 + status13

        # Column 1: Total, Water Column, Sedimented
        ax1 = axs[i, 0]
        ax1.fill_betweenx(y_status, water_col, color='blue', alpha=0.6, label='Water Column')
        ax1.fill_betweenx(y_status, sedimented, color='red', alpha=0.6, label='Sedimented')
        ax1.plot(total, y, color='k', linestyle='-', linewidth=2, label='Total')
        ax1.set_xlim(left=0)
        ax1.set_ylim(0, 430)
        ax1.invert_yaxis()
        ax1.grid(linestyle='--')
        if i == 0:
            ax1.set_title('Total Amount')
        ax1.set_ylabel(rf"$\bf{{{site}}}$" + "\nAvg. Depth [m]")
        if i == 8:
            ax1.set_xlabel("Particle Count")
        if i == 0:
            ax1.legend(fontsize=9, loc='lower right')

        # Column 2: Water Column Statuses
        ax2 = axs[i, 1]
        ax2.fill_betweenx(y_status, status1, color='black', alpha=0.6, label='Sewage')
        ax2.fill_betweenx(y_status, status2, color='red', alpha=0.6, label='Colloidal')
        ax2.fill_betweenx(y_status, status3, color='green', alpha=0.6, label='Marine')
        ax2.plot(water_col, y_status, color='k', linestyle='-', linewidth=2, label='Water Column')
        ax2.set_xlim(left=0)
        ax2.set_ylim(0, 430)
        ax2.invert_yaxis()
        ax2.grid(linestyle='--')
        if i == 0:
            ax2.set_title('Water Column')
            ax2.legend(fontsize=9, loc='lower right')
        if i == 8:
            ax2.set_xlabel("Particle  Count")

        # Column 3: Sedimented Statuses
        ax3 = axs[i, 2]
        ax3.fill_betweenx(y_status, status11, color='black', alpha=0.6, label='S. Sewage')
        ax3.fill_betweenx(y_status, status12, color='red', alpha=0.6, label='S. Colloidal')
        ax3.fill_betweenx(y_status, status13, color='green', alpha=0.6, label='S. Marine')
        ax3.plot(sedimented, y_status, color='k', linestyle='-', linewidth=2, label='Sedimented')
        ax3.set_xlim(left=0)
        ax3.set_ylim(0, 430)
        ax3.invert_yaxis()
        ax3.grid(linestyle='--')
        if i == 0:
            ax3.set_title('Sediment')
            ax3.legend(fontsize=9, loc='lower right')
        if i == 8:
            ax3.set_xlabel("Particle Count")

    fig.tight_layout()
    plt.show()

################ Status per region exlusive ############3
def status_states_regions_map(polygon_dict):
    water_statuses = [1, 2, 3]
    sediment_statuses = [11, 12, 13]

    # Dictionaries to hold separated data
    water_column = {}
    sedimented = {}

    for region_name, (_, df) in polygon_dict.items():
        valid_rows = df.dropna(subset=['lon', 'lat', 'depth', 'status'])

        if len(valid_rows) == 0:
            water_column[region_name] = {'lon': np.array([]), 'lat': np.array([]), 'depth': np.array([])}
            sedimented[region_name] = {'lon': np.array([]), 'lat': np.array([]), 'depth': np.array([])}
            continue

        # Concatenate all particle data
        lons = np.concatenate(valid_rows['lon'].values)
        lats = np.concatenate(valid_rows['lat'].values)
        depths = np.concatenate(valid_rows['depth'].values)
        statuses = np.concatenate(valid_rows['status'].values)

        # Apply status filters
        water_mask = np.isin(statuses, water_statuses)
        sed_mask = np.isin(statuses, sediment_statuses)

        # Store separated particle info
        water_column[region_name] = {
            'lon': lons[water_mask],
            'lat': lats[water_mask],
            'depth': depths[water_mask]
        }
        sedimented[region_name] = {
            'lon': lons[sed_mask],
            'lat': lats[sed_mask],
            'depth': depths[sed_mask]
        }
    return water_column, sedimented   
#
################# Calculation of volume in gridX and gridY ###########################
polygon_coords_N1 = [
    (Nx1_1[0], Ny1_1[0]),
    (Nx1_1[1], Ny1_1[1]),
    (Nx1_2[1], Ny1_2[1]),
    (Nx1_1[0], Ny1_2[0]),
    (Nx1_1[0], Ny1_1[0])
]
polygon_N1 = Polygon(polygon_coords_N1)
#Subregion Right
polygon_coords_N2 = [
    (Nx2_3[0], Ny2_3[0]),
    (Nx2_1[0], Ny2_1[0]),
    (Nx2_2[1], Ny2_2[1]),
    (Nx2_2[0], Ny2_2[0])
]
polygon_N2 = Polygon(polygon_coords_N2)
#Central Strait
polygon_coords_C1 = [
    (Cx1_1[1], Cy1_1[1]),
    (Cx1_1[0], Cy1_1[0]),
    (Cx1_2[0], Cy1_2[0]),
    (Cx1_2[1], Cy1_2[1]),
    (Cx1_3[1], Cy1_3[1]),
    (Cx1_4[1], Cy1_4[1])
]
polygon_C1 = Polygon(polygon_coords_C1)
#Southern Strait and Subregions
#Subregion North
polygon_coords_S1 = [
    (Sx1_1[0], Sy1_1[0]),
    (Sx1_1[1], Sy1_1[1]),
    (Sx1_2[1], Sy1_2[1]),
    (Sx1_5[1], Sy1_5[1]),
    (Sx1_6[0], Sy1_6[0]),
    (Sx1_4[0], Sy1_4[0])
]
polygon_S1 = Polygon(polygon_coords_S1)
#Subregion South
polygon_coords_S2 = [
    (Sx2_1[1], Sy2_1[1]),
    (Sx2_1[0], Sy2_1[0]),
    (Sx2_2[1], Sy2_2[1]),
    (Sx2_3[1], Sy2_3[1]),
    (Sx2_4[1], Sy2_4[1])
]
polygon_S2 = Polygon(polygon_coords_S2)
#Haro Strait
polygon_coords_H1 = [
    (Hx1_1[1], Hy1_1[1]),
    (Hx1_1[0], Hy1_1[0]),
    (Hx1_2[0], Hy1_2[0]),
    (Hx1_3[0], Hy1_3[0])

]
polygon_H1 = Polygon(polygon_coords_H1)
#Juan de Fuca Strait
polygon_coords_J1 = [
    (Jx1_1[0], Jy1_1[0]),
    (Jx1_1[1], Jy1_1[1]),
    (Jx1_2[1], Jy1_2[1]),
    (Jx1_3[1], Jy1_3[1])

]
polygon_J1 = Polygon(polygon_coords_J1)
#
# Fraser River
polygon_coords_F1 = [
    (Fx1_1[0], Fy1_1[0]),
    (Fx1_1[1], Fy1_1[1]),
    (Fx1_2[1], Fy1_2[1]),
    (Fx1_3[1], Fy1_3[1])
]
polygon_F1 = Polygon(polygon_coords_F1)
#**Howe Sound**
polygon_coords_HW1 = [
    (HWx1_1[1], HWy1_1[1]),
    (HWx1_1[0], HWy1_1[0]),
    (HWx1_2[1], HWy1_2[1]),
    (HWx1_3[1], HWy1_3[1]),
    (HWx1_4[1], HWy1_4[1]),
    (HWx1_5[1], HWy1_5[1]),
    
]
polygon_HW1 = Polygon(polygon_coords_HW1)
#
def volume_region(polygon):
    volume = xr.open_dataset('/ocean/vvalenzuela/MOAD/grid2/mesh_mask202108_TDV.nc')['volume']
    mask = xr.open_dataset(path['mask'])['tmask'][0]
    #
    x = volume['x']  
    y = volume['y'] 
    #
    xx, yy = np.meshgrid(x, y)
    #
    polygon_mask_2d = np.array([
        [polygon.contains(Point(xx[j, i], yy[j, i])) for i in range(len(x))]
        for j in range(len(y))
    ])

    #
    nz = volume.sizes['z']
    polygon_mask_3d = np.repeat(polygon_mask_2d[np.newaxis, :, :], nz, axis=0)

    #
    combined_mask = (polygon_mask_3d & (mask.values == 1))

    #
    mask_da = xr.DataArray(combined_mask, dims=volume.dims, coords=volume.coords)

    #
    volume_in_polygon_water = volume.where(mask_da).sum().item()
    #
    return volume_in_polygon_water
#
def volumes():
    regions_volumes = [volume_region(polygon_N1), volume_region(polygon_N2),
                       volume_region(polygon_C1), volume_region(polygon_S1), 
                       volume_region(polygon_S2), volume_region(polygon_H1),
                       volume_region(polygon_J1), volume_region(polygon_F1),
                       volume_region(polygon_HW1)]
    regions_string = ['N1', 'N2', 'C1', 'S1', 'S2', 'H1', 'J1', 'F1', 'HW1']
    return regions_volumes, regions_string