#!/usr/bin/env python3

"""This is a work in progress, its messy and not fully functional!"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters
from scipy import stats
import seaborn as sns #All sns features can be deleted without affecting program function



class schoolfield:
    def __init__(self, data, index):
        self.k = 8.62e-5 #Boltzmann constant
        self.Tref = 273.15 #Reference temperature - 0C
        self.data = self.clean_dataset(data)
        self.index = index #ID for the model
        
        self.temps = np.array(self.data['K'].values) #Temperatures
        self.responses = np.array(self.data['Cor_Trait_Value'].values) #Growth rates
        
        self.set_name()
        self.get_T_pk() #Assign value for T - peak
        self.estimate_E_init()
        self.estimate_B0()
        self.set_parameters()
        self.fit_model()
        self.smooth()
        self.get_AIC()
        
            
    def clean_dataset(self, data):
        "Normalise each dataset"
        #Transform temps to kelvin
        data['K'] = data['ConTemp'] + 273.15
        
        # Convert corrected value from s^-1 to d^-1
        data['Cor_Trait_Value'] = data['StandardisedTraitValue'] * 60 * 60 * 24 # Convert corrected value from s^-1 to d^-1
        
        #If any trait values are negative then subtract the smallest value to normalise
        minimum_trait_value  = data['Cor_Trait_Value'].min()
        
        if minimum_trait_value < 0:
            data['Cor_Trait_Value'] -= minimum_trait_value - 10E-10 #Get rid of any 0s
        else:
            data['Cor_Trait_Value'] += 10E-10
            
        return data
    
    def get_T_pk(self):
        "Find the temperature at which maximum growth is observed"
        self.Tpk_row = dataset['Cor_Trait_Value'].idxmax()
        self.T_pk = dataset.loc[self.Tpk_row]['K'] #Temperature at which rate is maximum
        
    def estimate_E_init(self):
        "Estimate energy value using the slope of the values to the peak"
        upslope = self.data.loc[:self.Tpk_row] #slice data so only values less than TKP row are included
        if upslope.shape[0] > 1:
            temps = upslope['K']
            responses = upslope['Cor_Trait_Value']
            x = 1 / (self.k * temps)
            y = np.log(responses)
            
            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
            
            self.E_init = abs(slope)
        else:
            self.E_init = 0.6 #Default value
            
    def estimate_B0(self):
        "Returns the response at the tempetature closest to Tref"
        if self.temps.min() > self.Tref:
            self.B0 = np.log(self.responses.min())
        else:
            self.B0 = np.log(self.responses[self.temps <= self.Tref].max())
            
    def set_name(self):
        "Set species name to be applied to plot title"
        genus = self.data['ConGenus'][0]
        species = self.data['ConSpecies'][0]

        self.name = ' '.join([str(self.index), ':' ,genus, species])
        
    def set_parameters(self):
        "Create a parameters object using out guesses, these will then be fitted using least squares regression"
        self.parameters = Parameters()
        #                   Name,      Start,   Can_Vary, Lower, Upper
        self.parameters.add_many(('B0_start', self.B0, True, -np.inf, np.inf,  None),
                          ('E', self.E_init, True, 0, np.inf,  None),
                          ('E_D',self.E_init * 4, True, 0, np.inf,  None),
                          ('T_pk', self.T_pk, True, 273.15-50, 273.15+150,  None))

    def fit(self, params, temps):
        "Fit a schoolfield curve to a list of temperature values"
        parameter_vals = params.valuesdict()
        B0 = parameter_vals['B0_start']
        E = parameter_vals['E']
        E_D = parameter_vals['E_D']
        T_pk = parameter_vals['T_pk']   

        fit = B0 + np.log(np.exp((-E / self.k) * ((1 / temps) - (1 / self.Tref))) /\
                        (1 + (E/(E_D - E)) * np.exp(E_D / self.k * (1 / self.T_pk - 1 / temps)))
                        )
        return fit
        
    def schoolfield_fit(self, params, temps, responses):
        "Called by fit model only, generates residuals using test values"
                        
        residuals = np.exp(self.fit(params, self.temps)) - responses
        
        return residuals
        
    def fit_model(self):
        "Least squares regression to minimise fit"
        self.model = minimize(self.schoolfield_fit, 
                              self.parameters, 
                              args=(self.temps, self.responses),
                              method="leastsq")
                              
        self.R2 = 1 - np.var(self.model.residual) / np.var(self.responses)

    def smooth(self):
        "Pass an interpolated list of temperature values back through the curve function to generate a smooth curve"
        self.smooth_x = np.arange(self.temps.min() - 1, self.temps.max() + 1, 0.1)
        self.smooth_y = np.exp(self.fit(self.model.params, self.smooth_x))
        
    def get_AIC(self):
        k = self.model.nvarys
        n = self.model.ndata
        rss = self.model.chisqr
        
        self.AIC = np.log((2 * np.exp(1)) / n) + n + 2 + n * np.log(rss) + 2 * k

    def __str__(self):
        vars = [self.name, self.B0, self.E_init, self.T_pk]
        text = """
        {0[0]}
        
        B0 = {0[1]:.2f}
        E = {0[2]:.2f}
        T Peak = {0[3]:.2f}
        
        -----------------------------------
        """.format(vars)
        return text
        
    def plot(self):
        f = plt.figure()
        sns.set(style="ticks")
        ax = f.add_subplot(111)
        plt.plot(self.temps, self.responses, marker='o', linestyle='None')
        plt.plot(self.smooth_x, self.smooth_y, marker='None')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Response')
        plt.title(self.name, fontsize=14, fontweight='bold')
        plt.text(0.1, 0.9,'R2: {0:.2f}\nAIC: {1:.2f}'.format(self.R2, self.AIC), ha='center', va='center', transform=ax.transAxes)
        sns.despine()
        
        plt.savefig('../results/{}.png'.format(self.index), bbox_inches='tight')
       
        
        
        
def get_datasets(path):
    "Create a set of temperature response curve datasets from csv"
    data = pd.read_csv(path, encoding = "ISO-8859-1") #Open in latin 1   
    ids = pd.unique(data['OriginalID']).tolist() #Get unique identifiers
    #create a dictionary of datasets for easy access later
    #Remove negative values if there are any to allow log transpforms
    Datasets = {}
    for id in ids:
        curve_data = data.loc[data['OriginalID'] == id] #seperate data by uniqueID
        curve_data = curve_data.sort_values(['ConTemp']).reset_index() #sort so rows are in temperature order, reset index to 0  
        Datasets[id] = curve_data
    return Datasets    
    
data_path = '../Data/Tom_Smith_IDs.csv'
Datasets = get_datasets(data_path)

for i in Datasets.keys():
    dataset = Datasets[i]
    if dataset.shape[0] > 3: #Must have more datapoints than number of variables
        model = schoolfield(dataset, i)
        print(model)
        model.plot()