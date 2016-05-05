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

data_path = '../Data/database.csv'
Datasets = get_datasets(data_path)
all_models = []

aux_parameters = ['FinalID', 'OriginalID', 'Citation', 'Latitude', 'Longitude', 'ConKingdom', 'ConPhylum', 'ConClass',
                  'ConOrder', 'ConFamily', 'ConGenus', 'ConSpecies', 'OptimalConditions', 'RespirationType']                
                  
for i in Datasets.keys():
    dataset = Datasets[i]
    if dataset.shape[0] > 3: #Must have more datapoints than number of variables
        est_params = estimate_parameters(dataset, i, aux_parameters)
        models = [schoolfield_two_factor(est_params, i), schoolfield_original_simple(est_params, i)]
        if dataset.shape[0] > 5: #This model has two additional variables
            models.append(schoolfield_original(est_params, i))
        
        all_models.append(models)
        best_model = max(models)
        if best_model:
            print(best_model)
            best_model.plot('../Results/fits')
         
all_models = rank_and_flatten(all_models)
output_csv(all_models, path = '../Results/summary.csv', aux_cols = aux_parameters)
output_csv(all_models, path = '../Data/summaries/summary.csv', aux_cols = aux_parameters)
    
print('Completed in: ', datetime.now() - starttime)