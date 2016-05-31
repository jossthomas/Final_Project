#!/usr/bin/env python3

"""
This is a very dirty script to find where the maximum response temperatures of archaea and bacteria actually lie.
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

aux_parameters = []       
model_names = ['Boltzmann_Arrhenius'] #'LM', 

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

Datasets = split_datasets(data, sep = 'ConKingdom', _sort = ['Est.Tpk'])  
mtg_bac = Datasets[0]
mtg_ach = Datasets[1]

def mesophile_thermophile(entry, temp):
    if entry > (temp + 273.15):
        return 'Thermophile'
    else:
        return 'Mesophile'

#create a mesophile_thermophile col
bac_max_e = 0
bac_max_e_tpk = 0
arch_max_e = 0
arch_max_e_tpk = 0

count = 1   
for temp in np.arange(20, 100, 0.5):
    bridge = lambda x: mesophile_thermophile(x, temp) #bridging function so I can send an argument to mesophile_thermophile
    
    mtg_bac['TempPref'] =  mtg_bac['Est.Tpk'].apply(bridge)
    mtg_ach['TempPref'] =  mtg_ach['Est.Tpk'].apply(bridge)

    TempPref_ach = split_datasets(mtg_ach, sep = 'TempPref', _sort = ['Est.Tpk'])
    TempPref_bac = split_datasets(mtg_bac, sep = 'TempPref', _sort = ['Est.Tpk'])

    all_temps = [('{}_archaea'.format(i), TempPref_ach[i]) for i in TempPref_ach.keys()] + [('{}_bacteria'.format(i), TempPref_bac[i]) for i in TempPref_bac.keys()]
    
                       
                       
    for group in all_temps:
        data = group[1]
        key = group[0]
        if data.shape[0] >= 20: 
            est_params = estimate_parameters(data, flags = param_est_flags, aux_parameters_names=['ConKingdom']) 
        
            model = fit_models(model_names, est_params, tag = count)[0]
            
            king = model.aux_parameters_values[0]
            mesothermo = model.species_name
            
            if king == 'Bacteria' and mesothermo == 'Thermophile':
                if model.final_E > bac_max_e:
                    bac_max_e = model.final_E
                    bac_max_e_tpk = temp
                    
            if king == 'Archaea' and mesothermo == 'Thermophile':
                if model.final_E > bac_max_e:
                    arch_max_e = model.final_E
                    arch_max_e_tpk = temp
            
        count += 1

print(bac_max_e_tpk, bac_max_e)
print(arch_max_e_tpk, arch_max_e)

print('Completed in: ', datetime.now() - starttime)