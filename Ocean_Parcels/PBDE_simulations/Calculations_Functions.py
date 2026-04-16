import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
#
def proportions_from_filename(filename):
    #
    data = xr.open_dataset(filename, engine = 'zarr')
    time_index = pd.to_datetime(data.time[0, :].values)

    # Categories based on status codes
    status_categories = {
        'Sewage W.C.': 1,
        'Colloidal W.C.': 2,
        'Marine W.C.': 3,
        'Sewage S.': 11,
        'Colloidal S.': 12,
        'Marine S.': 13,
        'Out JdF': 7,
        'Out Js': 8    }



    # Initialize results dict
    results = {label: [] for label in ['Initial'] + list(status_categories.keys())}
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



    # Convert counts to proportions (%)
    proportions = {
        label: (np.array(counts) / np.array(total_particles)) * 100
        for label, counts in results.items()
    }

    return pd.DataFrame(proportions, index=time_index)

###############################################################################
### REGIONS INFO ###
def metrics_table(filename, title = 'Simulation Metrics', plot = True):
    import Regions_functions_V2
    polygon_dict = Regions_functions_V2.polygon_definition(filename)
    status_vertical_N1 = Regions_functions_V2.vertical_status_profiles(polygon_dict['N1'], 80)
    status_vertical_S1 = Regions_functions_V2.vertical_status_profiles(polygon_dict['S1'], 80)
    status_vertical_H1 = Regions_functions_V2.vertical_status_profiles(polygon_dict['H1'], 80)
    volumes_regions_depths = Regions_functions_V2.volume_by_depth_all_regions(polygon_dict)
#
    def interpolate_volume_profile(volume_profile, number_of_depths):
        # Get the existing depth coordinate values (z-levels mapped to actual depth)
        depth_vals = volume_profile['depth'].values
        volume_vals = volume_profile.values

        # Filter out NaNs
        valid = ~np.isnan(depth_vals) & ~np.isnan(volume_vals)
        depth_vals = depth_vals[valid]
        volume_vals = volume_vals[valid]

        # Define new depth bins (higher resolution)
        new_depths = np.linspace(depth_vals.min(), depth_vals.max(), number_of_depths)

        # Interpolate depths
        interp_volume_vals = np.interp(new_depths, depth_vals, volume_vals)
        #
        interpolated_volume = xr.DataArray(
            interp_volume_vals,
            coords={'depth': new_depths},
            dims='depth',
            name='volume'
        )

        return interpolated_volume, new_depths
    ##################
    depth_bins_regions = 81
    #####
    volume_N1, depths_N1 = interpolate_volume_profile(volumes_regions_depths['N1'], depth_bins_regions)
    volume_N2, depths_N2 = interpolate_volume_profile(volumes_regions_depths['N2'], depth_bins_regions)
    volume_S1, depths_S1 = interpolate_volume_profile(volumes_regions_depths['S1'], depth_bins_regions)
    volume_SP, depths_SP = interpolate_volume_profile(volumes_regions_depths['SP'], depth_bins_regions)
    volume_HW1, depths_HW1 = interpolate_volume_profile(volumes_regions_depths['HW1'], depth_bins_regions)
    volume_S2, depths_S2 = interpolate_volume_profile(volumes_regions_depths['S2'], depth_bins_regions)
    volume_H1, depths_H1 = interpolate_volume_profile(volumes_regions_depths['H1'], depth_bins_regions)
    #################
    status_vertical_N1 = Regions_functions_V2.vertical_status_profiles_V2(polygon_dict['N1'], depth_bin_edges=depths_N1)
    status_vertical_N2 = Regions_functions_V2.vertical_status_profiles_V2(polygon_dict['N2'], depth_bin_edges=depths_N2)
    status_vertical_S1 = Regions_functions_V2.vertical_status_profiles_V2(polygon_dict['S1'], depth_bin_edges=depths_S1)
    status_vertical_SP = Regions_functions_V2.vertical_status_profiles_V2(polygon_dict['SP'], depth_bin_edges=depths_SP)
    status_vertical_HW1 = Regions_functions_V2.vertical_status_profiles_V2(polygon_dict['HW1'], depth_bin_edges=depths_HW1)
    status_vertical_S2 = Regions_functions_V2.vertical_status_profiles_V2(polygon_dict['S2'], depth_bin_edges=depths_S2)
    status_vertical_H1 = Regions_functions_V2.vertical_status_profiles_V2(polygon_dict['H1'], depth_bin_edges=depths_H1)
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
    def vertical_ratios(status_vertical):
        numerator = status_vertical['Particles Status 3'].astype(float)
        denominator = (
            status_vertical['Particles Status 3'] +
            status_vertical['Particles Status 2'] +
            status_vertical['Particles Status 1']
        ).astype(float)
        
        mask = (numerator == 0) | (denominator == 0)
        
        particulate_ratio = (numerator / denominator).mask(mask) * 100
        dissolved_ratio = 100 - particulate_ratio

        df = pd.DataFrame({
            'Depth': status_vertical['Avg. Depth'].values,
            'Dissolved (%)': dissolved_ratio.values,
            'Particulate (%)': particulate_ratio.values
        })

        return df
    N1_model, S1_model, SP_model = [
    vertical_ratios(status_vertical_N1), vertical_ratios(status_vertical_S1), vertical_ratios(status_vertical_SP)]
    #
    median_N1 = (N1_model['Particulate (%)'] / N1_model['Dissolved (%)']).median()
    median_SP = (SP_model['Particulate (%)'] / SP_model['Dissolved (%)']).median()
    median_S1 = (SP_model['Particulate (%)'] / SP_model['Dissolved (%)']).median()
    #
    model_ratios = [median_N1, median_SP, median_S1]
    #
    ratio_MC_Water = np.nanmedian(model_ratios)
    #
    ### RATIO M/C WATER ###
    #
    volumes_regions, _ = Regions_functions_V2.volumes()
    NSoG_OP = status_vertical_N2['Particles Status 11'].sum() + status_vertical_N2['Particles Status 12'].sum() + status_vertical_N2['Particles Status 13'].sum() 
    SSoGa_OP = status_vertical_S1['Particles Status 11'].sum() + status_vertical_S1['Particles Status 12'].sum() + status_vertical_S1['Particles Status 13'].sum()
    SSoGb_OP = status_vertical_S2['Particles Status 11'].sum() + status_vertical_S2['Particles Status 12'].sum() + status_vertical_S2['Particles Status 13'].sum() 
    IPlume_OP = status_vertical_SP['Particles Status 11'].sum() + status_vertical_SP['Particles Status 12'].sum() + status_vertical_SP['Particles Status 13'].sum() 
    HWSound_OP = status_vertical_HW1['Particles Status 11'].sum() + status_vertical_HW1['Particles Status 12'].sum() + status_vertical_HW1['Particles Status 13'].sum()  
    OP_model = [NSoG_OP / volumes_regions[1], SSoGa_OP / volumes_regions[4], IPlume_OP / volumes_regions[5], HWSound_OP / volumes_regions[6], SSoGb_OP / volumes_regions[8]]
    OP_modeli = OP_model / np.nansum(OP_model)    
    Concentration_Sediment = np.nanmedian(OP_modeli)
    #
    ### Simulation Info ###
    #
    data = xr.open_dataset(filename, engine='zarr')
    data_water = data.where((data.status > 0) & (data.status < 4))
    data_sediment = data.where(data.status > 10)
    #
    depth_mean_water = data_water.z.mean().values
    depth_mean_sediment = data_sediment.z.mean().values
    ####
    ####
    #if plot == True:
    #### METRICS TABLE ####
    #plt.rcParams.update({'font.size': 30})
        #column_labels = ['Simulation Label','M/C W. C.', 'M/C S.', 'N1/S1 Colloidal W. C.', 'H1/S1 Colloidal W. C.', 'Out JdF', 'Out Mixing Region', 'Colloidal W. C.', 'Mean Depth W. C.', 'Mean Depth S.']
        #
        #metrics = [title,
        #    f"{ratio_MC_Water:.5f}",
        #    f"{Concentration_Sediment:.5f}",
        #    f"{ratio_N1_S1:.5f}",
        #    f"{ratio_H1_S1:.5f}",       
        #    f"{proportions['Out JdF'].values[-1]:.5f} %",
        #    f"{proportions['Out Haro Mix'].values[-1]:.5f} %",
        #    f"{proportions['Colloidal Water'].values[-1]:.5f} %",
        #    f"{depth_mean_water:.5f} m",
        #    f"{depth_mean_sediment:.5f} m"
        #]
        #
        #fig, ax = plt.subplots(figsize=(16, 3))  
        #ax.axis("off")
        #
        #
        #table = ax.table(
        #    cellText=[metrics],
        #    colLabels=column_labels,
        #    loc="center",
        #    cellLoc="center"
        #)
        #
        #for (row, col), cell in table.get_celld().items():
        #    cell.set_height(0.3)
        #    cell.get_text().set_fontsize(50)
        #    if row == 0:
        #        cell.get_text().set_weight('bold')
        #
        #
        #table.scale(1.5, 3)
        #plt.tight_layout()
        #plt.show()    
    #
    df = pd.DataFrame([{
        'Simulation Label': title,
        'M/C': ratio_MC_Water,
        'Conc. S.': Concentration_Sediment,
        'N1/S1': ratio_N1_S1,
        'H1/S1': ratio_H1_S1,
        'Out JdF (%)': proportions['Out JdF'].values[-1],
        'Out Js (%)': proportions['Out Js'].values[-1],
        'Total Colloidal W. C. (%)': proportions['Colloidal W.C.'].values[-1],
        'W. C.': depth_mean_water,
        'Sed.': depth_mean_sediment
    }])    
    return df
    




















