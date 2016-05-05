#!/usr/bin/env python3

"""
This program provides a framework for and implementation of least squares fitting of various 
thermal response models based on experimental data

Written in Python 3.5 Anaconda Distribution
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, fit_report
from scipy import stats
import seaborn as sns #All sns features can be deleted without affecting program function
from General_Models import physiological_growth_model, Boltzmann_Arrhenius, schoolfield_two_factor, schoolfield_original_simple, schoolfield_original          
            
def get_datasets(path):
    "Create a set of temperature response curve datasets from csv"
    data = pd.read_csv(path, encoding = "ISO-8859-1") #Open in latin 1
    data['FinalID'] = pd.factorize(data['OriginalID'])[0]
    ids = pd.unique(data['FinalID']).tolist() #Get unique identifiers
    #create a dictionary of datasets for easy access later
    Datasets = {}
    for id in ids:
        curve_data = data.loc[data['FinalID'] == id] #seperate data by uniqueID
        curve_data = curve_data.sort_values(['ConTemp']).reset_index() #sort so rows are in temperature order, reset index to 0  
        Datasets[id] = curve_data
    return Datasets    
    
data_path = '../Data/methanogens.csv'
Datasets = get_datasets(data_path)
all_models = []

for i in Datasets.keys():
    print(i)
    dataset = Datasets[i]
    if dataset.shape[0] > 3: #Must have more datapoints than number of variables
        models = [schoolfield_two_factor(dataset, i), schoolfield_original_simple(dataset, i)]
        if dataset.shape[0] > 5: #This model has two additional variables
            models.append(schoolfield_original(dataset, i))
        
        best_model = max(models)
        all_models.append(best_model)

print('Models fitted')

#Create a blank dataframe
output = pd.DataFrame(columns=("Species", "Model_name", "Trait", "Latitude", "Longitude", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "B0", "E",
                               "E_D", "E_D_L", "Est.Tpk", "Est.Tmin", "Est.Tmax", "Max response", "R_Squared", "AIC", "BIC", "Temperature", "Response", "Observed")) 
#Add results to dataframe
print('output made')

iter = 0 
for model in all_models:
    print(model.name)
    smooth_x = list(model.smooth_x)
    smooth_y = list(model.smooth_y)
    x_obs = list(model.temps)
    y_obs = list(model.responses)
    final_E_D = getattr(model, 'final_E_D', "NA")
    final_E_D_L = getattr(model, 'final_E_D_L', "NA")
    final_T_H = getattr(model, 'final_T_H', "NA")
    final_T_H_L = getattr(model, 'final_T_H_L', "NA")
    for i, x in enumerate(smooth_x):
        temp = x
        response = smooth_y[i]
        
        output.loc[iter] = [model.name, model.model_name, model.trait, model.latitude, model.longditude, model.kingdom, model.phylum, model.class_, model.order, model.family, 
                            model.genus, model.final_B0, model.final_E, final_E_D, final_E_D_L, model.tpk_est, model.lower_percentile, model.upper_percentile, 
                            model.max_response_est, model.R2, model.AIC, model.BIC, temp, response, False]
        iter += 1
        
    for i, x in enumerate(x_obs):
        temp = x
        response = y_obs[i]
        
        output.loc[iter] = [model.name, model.model_name, model.trait, model.latitude, model.longditude, model.kingdom, model.phylum, model.class_, model.order, model.family, 
                            model.genus, model.final_B0, model.final_E, final_E_D, final_E_D_L, model.tpk_est, model.lower_percentile, model.upper_percentile, 
                            model.max_response_est, model.R2, model.AIC, model.BIC, temp, response, True]
        iter += 1
        
output = output.sort_values(['Species']).reset_index()
output.to_csv('../Data/summaries/methanogens_summary.csv')