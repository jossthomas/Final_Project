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
                           read_database, \
                           fit_models, \
                           bootstrap_model, \
                           split_datasets, \
                           rank_and_flatten, \
                           compile_models
import numpy as np
import pandas as pd
from datetime import datetime

starttime = datetime.now()

data_path = '../Data/summaries/aggregate_data.csv'
data = read_database(data_path) 

aux_parameters = [] #no supplementary information needed
analysis_levels = ['ConKingdom', 'ConPhylum', 'ConClass', 'ConOrder', 'ConFamily', 'ConGenus'] #levels at which models should be fitted
model_names = ['Boltzmann_Arrhenius', 'LM'] #models to fit

for level in analysis_levels:
    #set the axis type of the plot based on analysis level.
    if level == 'ConKingdom':
        hist_axes = True
    else:
        hist_axes = False
    
    #see readme for more details on this dictionary
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
    
    Datasets = split_datasets(data, sep = level, _sort = ['Est.Tpk'])      
    
    count = 1
    
    for i in Datasets.keys():
        dataset = Datasets[i]
        tmax_check = len(dataset['Est.Tpk'].dropna())
        if tmax_check >= 5: #Must have more datapoints than number of variables
            est_params = estimate_parameters(dataset, aux_parameters, flags = param_est_flags) 
            print(est_params)
            
            models = fit_models(model_names, est_params, tag = (i + 1))
            models = [bootstrap_model(model, est_params) for model in models] #boostrap each model
            all_models.append(models)
            
            plot_path = '../Results/Maxima_fits/{}/standard'.format(level)
            plot_path_log = '../Results/Maxima_fits/{}/log'.format(level)
            plot_path_arrh = '../Results/Maxima_fits/{}/arrh'.format(level)
            
            for i, model in enumerate(models):
                model.plot(plot_path, plot_residuals=True, hist_axes = hist_axes, fit_stats = False, convert_kelvin = True)
                model.plot(plot_path_log, scale_type = 'log', plot_residuals=True, hist_axes = hist_axes, fit_stats = False, convert_kelvin = False)
                model.plot(plot_path_arrh, scale_type = 'arrhenius', plot_residuals=True, hist_axes = hist_axes, fit_stats = False, convert_kelvin = True)
                
            count += 1

    all_models = rank_and_flatten(all_models)
    summary_path = '../Results/Maxima_fits/{}_summary.csv'.format(level)
    compile_models(all_models, path = summary_path, aux_cols = aux_parameters, sortby=['Species', 'Model_name'])

#------------------------ Imagine this is actually a new script, it may as well be ------------------------ 
#Analyse Metabolism within kingdoms - slightly gnarly, I should write something better

Datasets = split_datasets(data, sep = 'ConKingdom', _sort = ['Est.Tpk'])  
archaea = split_datasets(Datasets[1], sep = 'Best_Guess', _sort = ['Est.Tpk'])
bacteria = split_datasets(Datasets[0], sep = 'Best_Guess', _sort = ['Est.Tpk'])
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
    tmax_check = len(data['Est.Tpk'].dropna())
    
    if tmax_check >= 5: #Must have more datapoints than number of variables
        est_params = estimate_parameters(data, flags = param_est_flags, aux_parameters_names=['ConKingdom']) 
        print(est_params)
        
        models = fit_models(model_names, est_params, tag = (key))
        models = [bootstrap_model(model, est_params) for model in models]

        all_models.append(models)
        
        plot_path = '../Results/Maxima_fits/metabolism/standard'
        plot_path_log = '../Results/Maxima_fits/metabolism/log'
        plot_path_arrh = '../Results/Maxima_fits/metabolism/arrh'
        
        for i, model in enumerate(models):
            model.plot(plot_path, plot_residuals=False)
            model.plot(plot_path_log, scale_type = 'log', plot_residuals=False)
            model.plot(plot_path_arrh, scale_type = 'arrhenius', plot_residuals=False)
        
all_models = rank_and_flatten(all_models)
summary_path = '../Results/Maxima_fits/metabolism_summary.csv'
compile_models(all_models, path = summary_path, sortby=['Species', 'Model_name'], aux_cols=['ConKingdom'])

#------------------------ Imagine this is actually a new script again, it may as well be ------------------------ 
#Split by mesophile thermophile

mtg_bac = Datasets[0]
mtg_ach = Datasets[1]

def mesophile_thermophile(entry):
    if entry > 323.15:
        return 'Thermophile'
    else:
        return 'Mesophile'

#create a mesophile_thermophile col
mtg_bac['TempPref'] =  mtg_bac['Est.Tpk'].apply(mesophile_thermophile)
mtg_ach['TempPref'] =  mtg_ach['Est.Tpk'].apply(mesophile_thermophile)

TempPref_ach = split_datasets(mtg_ach, sep = 'TempPref', _sort = ['Est.Tpk'])
TempPref_bac = split_datasets(mtg_bac, sep = 'TempPref', _sort = ['Est.Tpk'])

all_temps = [('{}_archaea'.format(i), TempPref_ach[i]) for i in TempPref_ach.keys()] + [('{}_bacteria'.format(i), TempPref_bac[i]) for i in TempPref_bac.keys()]
all_models = []

param_est_flags = {'Trait_Col_Name': 'Trait',
                       'X_vals_col_name': 'Est.Tpk',
                       'Y_vals_col_name': 'Max.response',
                       'x_val_conversion': 1,
                       'log_data': False,
                       'is_celcius': False,
                       'species_data': False,
                       'full_strain_col_name': 'TempPref',
                       'genus_level_col_name': 'ConKingdom',
                       'species_level_col_name': 'TempPref'}

count = 1                      
                       
for group in all_temps:
    data = group[1]
    key = group[0]
    tmax_check = len(data['Est.Tpk'].dropna())
    
    if tmax_check >= 5: #Must have more datapoints than number of variables
        est_params = estimate_parameters(data, flags = param_est_flags, aux_parameters_names=['ConKingdom']) 
        print(est_params)
        
        models = fit_models(model_names, est_params, tag = count)
        models = [bootstrap_model(model, est_params) for model in models]
        all_models.append(models)
        
        plot_path = '../Results/Maxima_fits/TempPref/standard'
        plot_path_log = '../Results/Maxima_fits/TempPref/log'
        plot_path_arrh = '../Results/Maxima_fits/TempPref/arrh'
        
        for i, model in enumerate(models):
            model.plot(plot_path, plot_residuals=False)
            model.plot(plot_path_log, scale_type = 'log', plot_residuals=False)
            model.plot(plot_path_arrh, scale_type = 'arrhenius', plot_residuals=False)
            
        count += 1
        
all_models = rank_and_flatten(all_models)
summary_path = '../Results/Maxima_fits/Temp_Group_summary.csv'
compile_models(all_models, path = summary_path, sortby=['Species', 'Model_name'], aux_cols=['ConKingdom'])                      

print('Completed in: ', datetime.now() - starttime)