#!/usr/bin/env python3

"""
This program provides a framework for and implementation of least squares fitting of various 
thermal response models based on experimental data

Written in Python 3.5 Anaconda Distribution
"""
from Thermal_Models import estimate_parameters, \
                           physiological_growth_model, \
                           Boltzmann_Arrhenius, \
                           schoolfield_two_factor, \
                           schoolfield_original_simple, \
                           schoolfield_original, \
                           get_datasets, \
                           rank_and_flatten, \
                           output_csv
import numpy as np
import pandas as pd
from datetime import datetime

starttime = datetime.now()

data_path = '../Data/summaries/aggregate_data.csv'

aux_parameters = []       
analysis_levels = ['ConKingdom', 'ConPhylum', 'ConClass', 'ConOrder', 'ConFamily', 'ConGenus']

for level in analysis_levels:
    param_est_flags = {'Trait_Col_Name': 'Trait',
                       'X_vals_col_name': 'Est.Tpk',
                       'Y_vals_col_name': 'Max.response',
                       'x_val_conversion': 1,
                       'log_data': False,
                       'is_celcius': False,
                       'species_data': False,
                       'full_strain_col_name': level,
                       'genus_level_col_name': level,
                       'species_level_col_name': level}
    
    all_models = []
    
    Datasets = get_datasets(data_path, sep = level, _sort = ['Est.Tpk'])                            
    for i in Datasets.keys():
        dataset = Datasets[i]
        if dataset.shape[0] >= 6: #Must have more datapoints than number of variables
            est_params = estimate_parameters(dataset, aux_parameters, flags = param_est_flags) 
            print(est_params)
       
            models = [schoolfield_two_factor(est_params, i), Boltzmann_Arrhenius(est_params, i)]
            all_models.append(models)
            
            plot_path = '../Results/Maxima_fits/{}/standard'.format(level)
            plot_path_log = '../Results/Maxima_fits/{}/log'.format(level)
            plot_path_arrh = '../Results/Maxima_fits/{}/arrh'.format(level)
            
            for i, model in enumerate(models):
                model.plot(plot_path)
                model.plot(plot_path_log, scale_type = 'log')
                model.plot(plot_path_arrh, scale_type = 'arrhenius')

    all_models = rank_and_flatten(all_models)
    summary_path = '../Results/Maxima_fits/{}_summary.csv'.format(level)
    output_csv(all_models, path = summary_path, aux_cols = aux_parameters, sortby=['Species', 'Model_name'])
    
print('Completed in: ', datetime.now() - starttime)