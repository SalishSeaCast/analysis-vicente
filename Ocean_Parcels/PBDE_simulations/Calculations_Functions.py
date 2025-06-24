import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
#
def proportions_from_filename(filename):
    data = xr.open_dataset(filename)
    time_index = pd.to_datetime(data.time[0, :].values)
    
    # Fixed status codes
    categories = {
        'Sewage Water': 1,
        'Colloidal Water': 2,
        'Marine Water': 3,
        'Sewage Sediment': 11,
        'Colloidal Sediment': 12,
        'Marine Sediment': 13,
        'Out JdF': 7
    }

    #
    results = {label: [] for label in ['Initial'] + list(categories.keys())}
    total_particles = []

    for i in range(len(data.obs)):
        status_i = data.status[:, i].values
        valid = np.isfinite(status_i)

        total = np.count_nonzero(valid)
        total_particles.append(total)

        # Initial: any status < 0
        results['Initial'].append(np.count_nonzero(valid & (status_i < 0)))

        # Fixed-status categories
        for label, code in categories.items():
            results[label].append(np.count_nonzero(valid & (status_i == code)))

    # Calculate proportions
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
    ### Ratio N1-S1 ###
    # JUST COLLOIDAL
    sum_water_N1 = status_vertical_N1['Particles Status 2'].sum() #+ status_vertical_N1['Particles Status 3'].sum()
    sum_water_S1 = status_vertical_S1['Particles Status 2'].sum() #+ status_vertical_S1['Particles Status 3'].sum()
    ratio_N1_S1 = sum_water_N1 / sum_water_S1
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
    plt.rcParams.update({'font.size': 30})
    column_labels = ['Simulation Label','M/C Water Column', 'M/C Sediment', 'N1/S1 Colloidal Water Column', 'Status 7 (Out JdF)',
                'Total Colloidal Water Column', 'Mean Depth Water Column', 'Mean Depth Sediment']
    #
    metrics = [title,
        f"{ratio_MC_Water:.5f}",
        f"{ratio_MC_Sediment:.5f}",
        f"{ratio_N1_S1:.5f}",
        f"{proportions['Out JdF'].values[-1]:.5f} %",
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
        cell.get_text().set_fontsize(30)
        if row == 0:
            cell.get_text().set_weight('bold')
    
    #
    table.scale(1.2, 2.5)
    plt.tight_layout()
    plt.show()




















