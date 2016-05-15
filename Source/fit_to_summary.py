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
                           LM, \
                           get_datasets, \
                           rank_and_flatten, \
                           output_csv
import numpy as np
import pandas as pd
from datetime import datetime

starttime = datetime.now()

data_path = '../Data/summaries/aggregate_data.csv'
data = pd.read_csv(data_path, encoding = "ISO-8859-1") #Open in latin 1

aux_parameters = []       
analysis_levels = ['ConKingdom', 'ConPhylum', 'ConClass', 'ConOrder', 'ConFamily', 'ConGenus', 'Best_Guess']

                 
for level in analysis_levels:
    if level == 'ConKingdom':
        hist_axes = True
    else:
        hist_axes = False
        
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
    
    Datasets = get_datasets(data, sep = level, _sort = ['Est.Tpk'])      
    
    count = 1
    
    for i in Datasets.keys():
        dataset = Datasets[i]
        if dataset.shape[0] >= 5: #Must have more datapoints than number of variables
            est_params = estimate_parameters(dataset, aux_parameters, flags = param_est_flags) 
            print(est_params)
       
            models = [LM(est_params, count), Boltzmann_Arrhenius(est_params, count)]
            all_models.append(models)
            
            plot_path = '../Results/Maxima_fits/{}/standard'.format(level)
            plot_path_log = '../Results/Maxima_fits/{}/log'.format(level)
            plot_path_arrh = '../Results/Maxima_fits/{}/arrh'.format(level)
            
            for i, model in enumerate(models):
                model.plot(plot_path, plot_residuals=True, hist_axes = hist_axes)
                model.plot(plot_path_log, scale_type = 'log', plot_residuals=False, hist_axes = hist_axes)
                model.plot(plot_path_arrh, scale_type = 'arrhenius', plot_residuals=False, hist_axes = hist_axes)
                
            count += 1

    all_models = rank_and_flatten(all_models)
    summary_path = '../Results/Maxima_fits/{}_summary.csv'.format(level)
    output_csv(all_models, path = summary_path, aux_cols = aux_parameters, sortby=['Species', 'Model_name'])

#Analyse Metabolism within kingdoms - slightly knarly, I should write something better

Datasets = get_datasets(data, sep = 'ConKingdom', _sort = ['Est.Tpk'])  
archaea = get_datasets(Datasets[1], sep = 'Best_Guess', _sort = ['Est.Tpk'])
bacteria = get_datasets(Datasets[0], sep = 'Best_Guess', _sort = ['Est.Tpk'])
all_proks = [('{}_archaea'.format(i), archaea[i]) for i in archaea.keys()] + [('{}_bacteria'.format(i), bacteria[i]) for i in bacteria.keys()]
all_models = []

param_est_flags = {'Trait_Col_Name': 'Trait',
                       'X_vals_col_name': 'Est.Tpk',
                       'Y_vals_col_name': 'Max.response',
                       'x_val_conversion': 1,
                       'log_data': False,
                       'is_celcius': False,
                       'species_data': False,
                       'full_strain_col_name': 'Best_Guess',
                       'genus_level_col_name': 'ConKingdom',
                       'species_level_col_name': 'Best_Guess'}

for group in all_proks:
    data = group[1]
    key = group[0]
    if data.shape[0] >= 5: 
        est_params = estimate_parameters(data, flags = param_est_flags, aux_parameters_names=['ConKingdom']) 
        print(est_params)
        
        models = [LM(est_params, key), Boltzmann_Arrhenius(est_params, key)]
        all_models.append(models)
        
        plot_path = '../Results/Maxima_fits/metabolism/standard'
        plot_path_log = '../Results/Maxima_fits/metabolism/log'
        plot_path_arrh = '../Results/Maxima_fits/metabolism/arrh'
        
        for i, model in enumerate(models):
            model.plot(plot_path, plot_residuals=True)
            model.plot(plot_path_log, scale_type = 'log', plot_residuals=True)
            model.plot(plot_path_arrh, scale_type = 'arrhenius', plot_residuals=True)
        
all_models = rank_and_flatten(all_models)
summary_path = '../Results/Maxima_fits/metabolism_summary.csv'
output_csv(all_models, path = summary_path, sortby=['Species', 'Model_name'], aux_cols=['ConKingdom'])
    
print('Completed in: ', datetime.now() - starttime)