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
Datasets = split_datasets(data) #return a dictionary of growth curves split by ID

#names of the models I want to fit as strings
model_names = ['schoolfield_two_factor', 'schoolfield_original_simple', 'schoolfield_original']
all_models = [] #place to put fitted models

#Column names of parameters I want to use as explanatory variables when I analyse the data                  
aux_parameters = ['FinalID', 'OriginalID', 'Citation', 'Latitude', 'Longitude', 'ConKingdom', 'ConPhylum', 'ConClass',
                  'ConOrder', 'ConFamily', 'ConGenus', 'ConSpecies', 'OptimalConditions', 'Best_Guess'] 
                  
for i in Datasets.keys(): #Iterate through the dictionary by key
    dataset = Datasets[i] #get the growth curve
    est_params = estimate_parameters(dataset, aux_parameters) #Estimate starting parameters using regression
    models = fit_models(model_names, est_params, tag = i) #Fit 3 models in one line, returns a list containing them
    
    if models: #Check something has actually fitted
        best_model = max(models) #Find the best model
        if best_model.final_E > 0.001: #Checks that the model actually has a growth response
            all_models.append(models) #could also be all_models.append(best_model) if you aren't interested in the alternate fits
        print(best_model)
        best_model.plot('../Results/fits')
         
all_models = rank_and_flatten(all_models)
compile_models(all_models, path = '../Results/summary.csv', aux_cols = aux_parameters) #Save the summary in two places
compile_models(all_models, path = '../Data/summaries/summary.csv', aux_cols = aux_parameters)

print('Completed in: ', datetime.now() - starttime)