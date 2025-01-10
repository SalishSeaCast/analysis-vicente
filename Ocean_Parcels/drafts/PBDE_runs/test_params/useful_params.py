#### Useful parameters for PBDEs simulations ####
import xarray as xr
import numpy as np
import pandas as pd
#
#
### Percentage of each phase ###
#
#
def percentages(outfile_name):
    data = xr.open_zarr(outfile_name)
    #
    per_sewage_all = (len(np.where(data.status == 1.)[0])/(len(data.trajectory) * len(data.obs)))*100
    per_colloidal_all = (len(np.where(data.status == 2.)[0])/(len(data.trajectory) * len(data.obs)))*100    
    per_marine_all = (len(np.where(data.status == 3.)[0])/(len(data.trajectory) * len(data.obs)))*100
    per_bottom_all = (len(np.where(data.status == 4.)[0])/(len(data.trajectory) * len(data.obs)))*100
    #
    # INITIAL #
    per_sewage_initial = (len(np.where(data.status[:,0] == 1.)[0])/(len(data.trajectory)))*100
    per_colloidal_initial = (len(np.where(data.status[:,0] == 2.)[0])/(len(data.trajectory)))*100    
    per_marine_initial = (len(np.where(data.status[:,0] == 3.)[0])/(len(data.trajectory)))*100
    per_bottom_initial = (len(np.where(data.status[:,0] == 4.)[0])/(len(data.trajectory)))*100
    #
    # FINAL #
    per_sewage_final = (len(np.where(data.status[:,-1] == 1.)[0])/(len(data.trajectory)))*100
    per_colloidal_final = (len(np.where(data.status[:,-1] == 2.)[0])/(len(data.trajectory)))*100    
    per_marine_final = (len(np.where(data.status[:,-1] == 3.)[0])/(len(data.trajectory)))*100
    per_bottom_final = (len(np.where(data.status[:,-1] == 4.)[0])/(len(data.trajectory)))*100    
    #
    percentages_pbdes = np.array([per_sewage_all, per_colloidal_all, per_marine_all, per_bottom_all])
    percentages_initial = np.array([per_sewage_initial, per_colloidal_initial, per_marine_initial, per_bottom_initial])
    percentages_final = np.array([per_sewage_final, per_colloidal_final, per_marine_final, per_bottom_final])


    return percentages_initial, percentages_final, percentages_pbdes
#
#
def proportions(outfile_name):
    #
    data = xr.open_zarr(outfile_name)
    #
    colloidal = []
    marine = []
    bottom = []
    sewage = []
    #
    for i in range(len(data.obs)):
        len_1 = len(np.where(data.status[:,i] == 1.)[0])
        len_2 = len(np.where(data.status[:,i] == 2.)[0])
        len_3 = len(np.where(data.status[:,i] == 3.)[0])
        len_4 = len(np.where(data.status[:,i] == 4.)[0])
        #
        sewage.append(len_1)
        colloidal.append(len_2)
        marine.append(len_3)
        bottom.append(len_4)
    #    
    proportion_colloidal = np.array(colloidal)/len(data.trajectory)
    proportion_marine = np.array(marine)/len(data.trajectory)
    proportion_bottom = np.array(bottom)/len(data.trajectory)
    proportion_sewage = np.array(sewage)/len(data.trajectory)
    #
    # store data in pandas dataframe #
    pbdes_proportions = pd.DataFrame(np.transpose((proportion_sewage, proportion_colloidal, proportion_marine, proportion_bottom)), columns = ['Sewage', 'Colloidal', 'Marine', 'Bottom'], index = data.time[0,:].values)
    #
    return pbdes_proportions

