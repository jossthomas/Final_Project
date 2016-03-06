#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import Minimizer, minimize, Parameters, Parameter, report_fit, fit_report
from scipy import stats

class schoolfield:
    pass

def create_curve(temps, parameters):
    "Use the final parameters to create a smooth curve for plotting"
    temps = np.array([min(temps):max(temps):0.1])
    parameter_vals = params.valuesdict()
    B0 = parameter_vals['B0_start']
    E = parameter_vals['E']
    E_D = parameter_vals['E_D']
    T_pk = parameter_vals['T_pk']
    
    return(schoolfield_equation(temps, B0, E, E_D, T_pk))
    
    
def schoolfield_equation(temps, B0, E, E_D, T_pk, Tref=273.15, k=8.62e-5 ):
    eq = (B0 + np.log(np.exp((-E/k) * ((1/temps)-(1/Tref)))/(1 + (E/(E_D - E)) * np.exp(E_D/k * (1/T_pk - 1/temps)))))
    return residuals
    
def schoolfield_fit(params, temps, responses):
    #Unpack the parameters
    parameter_vals = params.valuesdict()
    B0 = parameter_vals['B0_start']
    E = parameter_vals['E']
    E_D = parameter_vals['E_D']
    T_pk = parameter_vals['T_pk']
    
    #Calculate the residuals of the model
    residuals = schoolfield_equation(temps, B0, E, E_D, T_pk)
    residuals = np.exp(residuals) - responses
    return residuals
    
def schoolfield_model(data, parameters):
    temps = np.array(data['K'].values)
    responses = np.array(data['OriginalTraitValue'].values)
    
    model = minimize(schoolfield_fit, parameters, args=(temps, responses),method="leastsq")
    
    output_parameters = model.params
    r_squared = 1-model.residual.var()/np.var(responses)
    n_vars = model.nvarys
    n_data = model.ndata
    chi_squared = model.chisqr
    aic = model.aic
    bic = model.bic
    model_vals = responses + model.residual
    
    return([model_vals, output_parameters, r_squared, n_vars, n_data, chi_squared, aic, bic])
    
def setup_parameters(B0_init, E_init, T_pk_init, theta=7):
    "Create a parameters object for the dataset which can then be fitted to the model"
    res = Parameters()
    res.add_many(('B0_start', B0_init, True, 0, 1000,  None),
                ('E', E_init, True, 0, 100,  None),
                ('E_D',4 *  E_init, True, 0, 100,  None),
                ('T_pk', T_pk_init, True, 273.15-50, 273.15+150,  None))
    return res
    
def estimate_E(data, T_PK_index, k=8.62e-5):
    "Estimate energy value using the slope of the values to the peak"
    upslope = data.loc[:T_PK_index]
    if upslope.shape[0] > 1:
        temps = upslope['K']
        responses = upslope['OriginalTraitValue']
        
        x = 1 / (k * temps)
        y = np.log(responses)
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
        return abs(slope)
    return 0.6 #default value
    
def estimate_B0(data, Tref=273.15):
    "Returns the rate at the tempetature closest to Tref"
    temps = data['K']
    responses = data['OriginalTraitValue']
    
    if min(temps) > Tref:
        return np.log(min(responses))
    else:
        rates_max_T = data.loc[temps <= Tref]['OriginalTraitValue']
        return np.log(max(rates_max_T))
    
def get_datasets(path):
    "Create a set of temperature response curve datasets from csv"
    data = pd.read_csv(path, encoding = "ISO-8859-1") #Open in latin 1   
    ids = pd.unique(data['OriginalID']).tolist() #Get unique identifiers
    data['K'] = data['ConTemp'] + 273.15 #Transform temps to kelvin
    
    #create a dictionary of datasets for easy access later
    #Remove negative values if there are any to allow log transpforms
    Datasets = {}
    for id in ids:
        curve_data = data.loc[data['OriginalID'] == id] #seperate data by uniqueID
        curve_data = curve_data.sort_values(['K']).reset_index() #sort so rows are in temperature order, reset index to 0  
        
        minimum_trait_value = min(curve_data['OriginalTraitValue']) #If minimum value < 0 then subtract from all values
        if minimum_trait_value < 0:
            curve_data['OriginalTraitValue'] -= minimum_trait_value - 10E-10 #Get rid of any 0s
        else:
            curve_data['OriginalTraitValue'] += 10E-10
            
        Datasets[id] = curve_data
    return Datasets

def plot_models(name, temps, responses, model):
    plt.figure()
    fig = plt.subplot(111)
    fig.plot(temps, responses, marker='o', color='b', linestyle='None', label='Result')
    fig.plot(temps, model, marker='None', color='r', label='Fit')
    plt.xlabel('Tempetature (K)')
    plt.ylabel('Response')
    plt.legend()
    plt.savefig('../results/{}.png'.format(name), bbox_inches='tight')

data_path = '../Data/Tom_Smith_IDs.csv'
Datasets = get_datasets(data_path)

for i in Datasets.keys():
    dataset = Datasets[i]
    if dataset.shape[0] > 3: #Must have more datapoints than number of variables
        print(i)
        temps = np.array(dataset['K'].values)
        responses = np.array(dataset['OriginalTraitValue'].values)
        
        T_PK_index = dataset['OriginalTraitValue'].idxmax()# Row of max trait value
        T_H_init = dataset.loc[T_PK_index]['K'] #Temperature at which rate is maximum
        E_init = estimate_E(dataset, T_PK_index) #Starting value for E
        B_init = estimate_B0(dataset) #Estimate of baseline metabolic rate
        Params = setup_parameters(B_init, E_init, T_H_init)
        result = schoolfield_model(dataset, Params)
        plot_models(i, temps, responses, result[0])


