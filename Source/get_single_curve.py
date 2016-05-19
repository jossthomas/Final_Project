from Thermal_Models import estimate_parameters, \
                           physiological_growth_model, \
                           Boltzmann_Arrhenius, \
                           schoolfield_two_factor, \
                           schoolfield_original_simple, \
                           schoolfield_original, \
                           read_database, \
                           fit_models, \
                           split_datasets, \
                           rank_and_flatten, \
                           compile_models
import numpy as np
import pandas as pd
from datetime import datetime

starttime = datetime.now()

data_path = '../Data/database.csv'
data = read_database(data_path) 

Datasets = split_datasets(data)
model_names = ['schoolfield_original']
all_models = []

for i in Datasets.keys():
    dataset = Datasets[i]
    est_params = estimate_parameters(dataset)
    if est_params.species_name == 'Halomonas subglaciescola':
        models = fit_models(model_names, est_params, tag = i)
        all_models.append(models)
         
all_models = rank_and_flatten(all_models)
compile_models(all_models, path = 'single.csv', whole_curves = True)

print('Completed in: ', datetime.now() - starttime)