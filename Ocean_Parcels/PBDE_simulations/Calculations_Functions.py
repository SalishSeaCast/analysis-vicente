import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
#
def proportions_from_filename(filename):
    import xarray as xr
    import numpy as np
    import pandas as pd

    data = xr.open_dataset(filename)
    time_index = pd.to_datetime(data.time[0, :].values)

    # Categories based on status codes
    status_categories = {
        'Sewage Water': 1,
        'Colloidal Water': 2,
        'Marine Water': 3,
        'Sewage Sediment': 11,
        'Colloidal Sediment': 12,
        'Marine Sediment': 13,
        'Out JdF': 7
    }

    # Extra categories using spatial filters
    spatial_categories = {
        'Haro Mix South': lambda lat, lon: (lat <= 48.46) & (lon <= -122) & (lon > -124.66),
        'Haro Mix North': lambda lat, lon: (lat < 48.7) & (lat > 48.46) & (lon > -124.66) & (lon < -124)
    }

    # Initialize results dict
    results = {label: [] for label in ['Initial'] + list(status_categories.keys()) + ['Out Haro Mix']}
    total_particles = []

    for i in range(len(data.obs)):
        status_i = data.status[:, i].values
        lat_i = data.lat[:, i].values
        lon_i = data.lon[:, i].values

        valid_status = np.isfinite(status_i)
        valid_pos = np.isfinite(lat_i) & np.isfinite(lon_i)

        total = np.count_nonzero(valid_status)
        total_particles.append(total)

        # Initial particles: status < 0
        results['Initial'].append(np.count_nonzero(valid_status & (status_i < 0)))

        # Status categories
        for label, code in status_categories.items():
            results[label].append(np.count_nonzero(valid_status & (status_i == code)))

        # Spatial categories
        mix_0 = np.count_nonzero(valid_pos & spatial_categories['Haro Mix South'](lat_i, lon_i))
        mix_1 = np.count_nonzero(valid_pos & spatial_categories['Haro Mix North'](lat_i, lon_i))
        
        results['Out Haro Mix'].append(mix_0 + mix_1)

    # Convert counts to proportions (%)
    proportions = {
        label: (np.array(counts) / np.array(total_particles)) * 100
        for label, counts in results.items()
    }

    return pd.DataFrame(proportions, index=time_index)

###############################################################################
### REGIONS INFO ###
def metrics_table(filename, title = 'Simulation Metrics'):
    import Regions_functions_V2
    polygon_dict = Regions_functions_V2.polygon_definition(filename)
    status_vertical_N1 = Regions_functions_V2.vertical_status_profiles(polygon_dict['N1'], 80)
    status_vertical_S1 = Regions_functions_V2.vertical_status_profiles(polygon_dict['S1'], 80)
    status_vertical_H1 = Regions_functions_V2.vertical_status_profiles(polygon_dict['H1'], 80)
    ### Ratio N1-S1 ###
    # JUST COLLOIDAL
    sum_water_N1 = status_vertical_N1['Particles Status 2'].sum() #+ status_vertical_N1['Particles Status 3'].sum()
    sum_water_S1 = status_vertical_S1['Particles Status 2'].sum() #+ status_vertical_S1['Particles Status 3'].sum()
    ratio_N1_S1 = sum_water_N1 / sum_water_S1
    ### Ratio H1-S1 ###
    # JUST COLLOIDAL
    sum_water_H1 = status_vertical_H1['Particles Status 2'].sum() #+ status_vertical_H1['Particles Status 3'].sum()
    ratio_H1_S1 = sum_water_H1 / sum_water_S1    
    ###############################################################################
    ### Proportions INFO ###
    #
    proportions = proportions_from_filename(filename)
    #
    ### RATIO M/C WATER ###
    #
    ratio_MC_Water = proportions['Marine Water'].values[-1] / proportions['Colloidal Water'].values[-1]
    #
    ### RATIO M/C WATER ###
    #
    ratio_MC_Sediment = proportions['Marine Sediment'].values[-1] / proportions['Colloidal Sediment'].values[-1]
    #
    ### Simulation Info ###
    #
    data = xr.open_dataset(filename)
    data_water = data.where((data.status > 0) & (data.status < 4))
    data_sediment = data.where(data.status > 10)
    #
    depth_mean_water = data_water.z.mean().values
    depth_mean_sediment = data_sediment.z.mean().values
    ####
    ####
    #### METRICS TABLE ####
    #plt.rcParams.update({'font.size': 30})
    column_labels = ['Simulation Label','M/C Water Column', 'M/C Sediment', 'N1/S1 Colloidal Water Column', 'H1/S1 Colloidal Water Column', 'Status 7 (Out JdF)', 'Out Mixing Region', 'Total Colloidal Water Column', 'Mean Depth Water Column', 'Mean Depth Sediment']
    #
    metrics = [title,
        f"{ratio_MC_Water:.5f}",
        f"{ratio_MC_Sediment:.5f}",
        f"{ratio_N1_S1:.5f}",
        f"{ratio_H1_S1:.5f}",       
        f"{proportions['Out JdF'].values[-1]:.5f} %",
        f"{proportions['Out Haro Mix'].values[-1]:.5f} %",
        f"{proportions['Colloidal Water'].values[-1]:.5f} %",
        f"{depth_mean_water:.5f} m",
        f"{depth_mean_sediment:.5f} m"
    ]
    #
    fig, ax = plt.subplots(figsize=(16, 3))  
    ax.axis("off")
    #
    #
    table = ax.table(
        cellText=[metrics],
        colLabels=column_labels,
        loc="center",
        cellLoc="center"
    )
    #
    for (row, col), cell in table.get_celld().items():
        cell.set_height(0.3)
        cell.get_text().set_fontsize(50)
        if row == 0:
            cell.get_text().set_weight('bold')
    
    #
    table.scale(1.5, 3)
    plt.tight_layout()
    plt.show()




















