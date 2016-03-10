#!/usr/bin/env python3

"""This is a work in progress, its messy and not fully functional!"""

import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, fit_report
from scipy import stats
import seaborn as sns #All sns features can be deleted without affecting program function

class growth_model:
    """
    Base class, not to be directly called. All other models will inherit from this!
    ALL child classes must have correctly named methods (fit, get_final_values)
    """
    
    k = 8.62e-5 #Boltzmann constant
    Tref = 273.15 #Reference temperature - 0C
    def __init__(self, data, index):
        self.data = self.clean_dataset(data)
        self.index = index #ID for the model
        
        self.temps = np.array(self.data['K'].values) #Temperatures
        self.responses = np.array(self.data['Cor_Trait_Value'].values) #Growth rates
        
        self.set_name()
        self.get_T_pk() #Assign value for T - peak
        self.estimate_T_H()
        self.estimate_T_H_L()
        self.estimate_E_init()
        self.estimate_B0()
        
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
        
    def estimate_T_H(self):
        "Find the temperature at which half the enzyme is inactivated by high temperatures"
        downslope = self.data.loc[self.Tpk_row:] #slice data so only the downwards slope is included (including TPK value)
        if downslope.shape[0] > 1:
            x = downslope['K']
            y = downslope['Cor_Trait_Value']
            
            #Not a linear thing, but this is fast
            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y) 
            x_intercept = -intercept/slope

            #find the value of K half way between T growth_max and T growth_0 (high), using 3 not 2 seems to work, no idea why
            self.T_H = self.T_pk + ((x_intercept - self.T_pk) / 3) 
        else:
            self.T_H = self.T_pk + 5 #Totally arbitary, probably wont fit

    def estimate_T_H_L(self):
        "Find the temperature at which half the enzyme is inactivated by high temperatures"
        upslope = self.data.loc[:self.Tpk_row] #slice data so only the downwards slope is included (including TPK value)
        if upslope.shape[0] > 1:
            x = upslope['K']
            y = upslope['Cor_Trait_Value']
            
            #Not a linear thing, but this is fast
            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y) 
            x_intercept = -intercept/slope

            #find the value of K half way between T growth_max and T growth_0 (low)
            self.T_H_L = self.T_pk - ((self.T_pk - x_intercept) / 2) 
        else:
            self.T_H_L = self.T_pk - 5 #Totally arbitary,  probably wont fit
            
    def estimate_E_init(self):
        "Estimate energy value using the slope of the values to the peak of an arrhenius plot"
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

        self.name = r' '.join([genus, species])
        
    def assess_model(self):
        k = self.model.nvarys
        n = self.model.ndata
        rss = sum(np.power(self.model.residual, 2))
        
        self.AIC = n * np.log((2 * np.pi) / n) + n + 2 + n * np.log(rss) + 2 * k
        self.BIC = n + n * np.log(2 * np.pi) + n * np.log(rss / n) + (np.log(n)) * (k + 1)
           
    def smooth(self):
        "Pass an interpolated list of temperature values back through the curve function to generate a smooth curve"
        self.smooth_x = np.arange(self.temps.min() - 3, self.temps.max() + 3, 0.1)
        self.smooth_y = np.exp(self.fit(self.model.params, self.smooth_x))
         
    def plot(self):
        textdata = [self.final_E,self.final_E_D, self.R2, self.AIC, self.BIC]
        title = '{}: {}'.format(self.index, self.name)
        
        f = plt.figure()
        sns.set_style("ticks", {'axes.grid': True})
        ax = f.add_subplot(111)
        
        plt.plot(self.smooth_x, self.smooth_y, marker='None', color='royalblue', linewidth=3)
        plt.plot(self.temps, self.responses, marker='o', linestyle='None', color='green')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Response')
        plt.title(title, fontsize=14, fontweight='bold')
        plt.text(0.05, 0.85,'E:  {0[0]:.2f}\nED: {0[1]:.2f}\nR2:  {0[2]:.2f}\nAIC: {0[3]:.2f}\nBIC: {0[4]:.2f}'.format(textdata),
                 ha='left', va='center', transform=ax.transAxes, color='darkslategrey')
        sns.despine()
        
        plt.savefig('../results/{}.png'.format(self.index), bbox_inches='tight')
        plt.close()


class schoolfield_two_factor(growth_model):
    def __init__(self, data, index):
        super().__init__(data, index) #Run the __init__ method from the base class
        self.index = str(index) + "_Sch_TF" #ID for the model
        self.set_parameters()
        self.fit_model()
        self.smooth()
        self.assess_model()
        self.get_final_values()
        
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
                        (1 + (E/(E_D - E)) * np.exp(E_D / self.k * (1 / T_pk - 1 / temps)))
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
        
    def get_final_values(self):
        "Get the final fitted values for the model"
        values = self.model.params.valuesdict()
        self.final_B0 = values['B0_start']
        self.final_E = values['E']
        self.final_E_D = values['E_D']
        self.final_T_pk = values['T_pk']   
        
    def __str__(self):
        "Allows print() to be called on the object"
        vars = [self.name, self.B0, self.final_B0, self.E_init, self.final_E, self.T_pk, self.final_T_pk,
                self.E_init * 4, self.final_E_D, self.R2, self.AIC, self.BIC]
        text = """\
        ---Schoolfield Two Factor Model---
        {0[0]}
        
        B0 est = {0[1]:.2f}
        B0 final = {0[2]:.2f}
        
        E est = {0[3]:.2f}
        E final = {0[4]:.2f}
        
        T Peak est = {0[5]:.2f}
        T Peak final =  {0[6]:.2f}
        
        E_D est = {0[7]:.2f}
        E_D final = {0[8]:.2f}
        
        R2: = {0[9]:.2f}
        AIC = {0[10]:.2f}
        BIC = {0[11]:.2f}
        
        -----------------------------------
        """.format(vars)
        return text

class schoolfield_original_simple(growth_model):
     # Original Sharpe-Schoolfield's model with low-temp inactivation term removed, very similar to previous model
     # This seems to generate pretty much identical results to the two factor model, but throws its toys out the pram sometimes.
     # In this version T_H is the temperature at which half the enzyme is denatured by heat stress
    def __init__(self, data, index):
        super().__init__(data, index) #Run the __init__ method from the base class
        self.index = str(index) + "_Sch_OS"
        self.set_parameters()
        self.fit_model()
        self.smooth()
        self.assess_model()
        self.get_final_values()
        
    def set_parameters(self):
        "Create a parameters object using out guesses, these will then be fitted using least squares regression"
        self.parameters = Parameters()
        #                   Name,      Start,   Can_Vary, Lower, Upper
        self.parameters.add_many(('B0_start', self.B0, True, -np.inf, np.inf,  None),
                          ('E', self.E_init, True, 0, np.inf,  None),
                          ('E_D',self.E_init * 4, True, 0, np.inf,  None),
                          ('T_H', self.T_H, True, self.T_pk, 273.15+170,  None))

    def fit(self, params, temps):
        "Fit a schoolfield curve to a list of temperature values"
        parameter_vals = params.valuesdict()
        B0 = parameter_vals['B0_start']
        E = parameter_vals['E']
        E_D = parameter_vals['E_D']
        T_H = parameter_vals['T_H']   
        
        fit = B0 + np.log(np.exp((-E / self.k) * ((1 / temps) - (1 / self.Tref)))\
                        /(1 + np.exp((E_D / self.k) * (1 / T_H - 1 / temps))))
        
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
        
    def get_final_values(self):
        "Get the final fitted values for the model"
        values = self.model.params.valuesdict()
        self.final_B0 = values['B0_start']
        self.final_E = values['E']
        self.final_E_D = values['E_D']
        self.final_T_H = values['T_H']   
        
    def __str__(self):
        "Allows print() to be called on the object"
        vars = [self.name, self.B0, self.final_B0, self.E_init, self.final_E, self.T_pk, self.T_H, 
                self.final_T_H, self.E_init * 4, self.final_E_D, self.R2, self.AIC, self.BIC]
        text = """\
        ---Schoolfield Original Model With Single T_H---
        {0[0]}
        
        B0 est = {0[1]:.2f}
        B0 final = {0[2]:.2f}
        
        E est = {0[3]:.2f}
        E final = {0[4]:.2f}
        
        TPK = {0[5]:.2f}
        T H est = {0[6]:.2f}
        T H final =  {0[7]:.2f}
        
        E_D est = {0[8]:.2f}
        E_D final = {0[9]:.2f}
        
        R2: = {0[10]:.2f}
        AIC = {0[11]:.2f}
        BIC = {0[12]:.2f}
     
        """.format(vars)
        return text        


class schoolfield_original(growth_model):
     # Original Sharpe-Schoolfield's model with low-temp inactivation term removed, very similar to previous model
     # In this version T_H_L is a low temperature enzyme inactivation constant (as if having a high temp one wasn't fun enough already)
    def __init__(self, data, index):
        super().__init__(data, index) #Run the __init__ method from the base class
        self.index = str(index) + "_Sch_O"
        self.set_parameters()
        self.fit_model()
        self.smooth()
        self.assess_model()
        self.get_final_values()
        
    def set_parameters(self):
        "Create a parameters object using out guesses, these will then be fitted using least squares regression, note additional T_H_L parameter"
        self.parameters = Parameters()
        #                   Name,      Start,   Can_Vary, Lower, Upper
        self.parameters.add_many(('B0_start', self.B0, True, -np.inf, np.inf,  None),
                          ('E', self.E_init, True, 0, np.inf,  None),
                          ('E_D',self.E_init * 4, True, 0, np.inf,  None),
                          ('E_D_L',self.E_init * (-6), True, -np.inf, 0,  None),
                          ('T_H', self.T_H, True, self.T_pk, 273.15+170,  None),
                          ('T_H_L', self.T_H_L, True, 273.15-70, self.T_pk,  None))

    def fit(self, params, temps):
        "Fit a schoolfield curve to a list of temperature values"
        parameter_vals = params.valuesdict()
        B0 = parameter_vals['B0_start']
        E = parameter_vals['E']
        E_D = parameter_vals['E_D']
        E_D_L = parameter_vals['E_D_L']
        T_H = parameter_vals['T_H']
        T_H_L = parameter_vals['T_H_L']         
        
        fit = B0 + np.log(np.exp((-E / self.k) * ((1 / temps) - (1 / self.Tref)))\
                        /(1 + np.exp((E_D_L / self.k) * (1 / T_H_L - 1 / temps)) + \
                              np.exp((E_D / self.k) * (1 / T_H - 1 / temps))))
        
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
        
    def get_final_values(self):
        "Get the final fitted values for the model"
        values = self.model.params.valuesdict()
        self.final_B0 = values['B0_start']
        self.final_E = values['E']
        self.final_E_D = values['E_D']
        self.final_E_D_L = values['E_D_L']
        self.final_T_H = values['T_H']
        self.final_T_H_L = values['T_H_L']           
        
    def __str__(self):
        "Allows print() to be called on the object"
        vars = [self.name, self.B0, self.final_B0, self.E_init, self.final_E, self.E_init * 4, self.final_E_D, self.E_init * (-6), self.final_E_D_L,
                self.T_pk, self.T_H, self.final_T_H, self.T_H_L, self.final_T_H_L, self.R2, self.AIC, self.BIC]
        text = """\
        ---Schoolfield Original Model With T_H and T_H_L---
        {0[0]}
        
        B0 est = {0[1]:.2f}
        B0 final = {0[2]:.2f}
        
        E est = {0[3]:.2f}
        E final = {0[4]:.2f}     
        
        E D est = {0[5]:.2f}
        E D final = {0[6]:.2f}
        
        E D L est = {0[7]:.2f}
        E D L final = {0[8]:.2f}
        
        TPK = {0[9]:.2f}
        T H est = {0[10]:.2f}
        T H final =  {0[11]:.2f}
        T H L est = {0[12]:.2f}
        T H L final =  {0[13]:.2f}  
        
        R2: = {0[14]:.2f}
        AIC = {0[15]:.2f}
        BIC = {0[16]:.2f}
     
        """.format(vars)
        return text                
        
def get_datasets(path):
    "Create a set of temperature response curve datasets from csv"
    data = pd.read_csv(path, encoding = "ISO-8859-1") #Open in latin 1   
    ids = pd.unique(data['OriginalID']).tolist() #Get unique identifiers
    #create a dictionary of datasets for easy access later
    Datasets = {}
    for id in ids[:15]:
        curve_data = data.loc[data['OriginalID'] == id] #seperate data by uniqueID
        curve_data = curve_data.sort_values(['ConTemp']).reset_index() #sort so rows are in temperature order, reset index to 0  
        Datasets[id] = curve_data
    return Datasets    
    
data_path = '../Data/Tom_Smith_IDs.csv'
Datasets = get_datasets(data_path)

for i in Datasets.keys():
    dataset = Datasets[i]
    if dataset.shape[0] > 3: #Must have more datapoints than number of variables
        model_TF = schoolfield_two_factor(dataset, i)
        print(model_TF)
        model_TF.plot()
        model_OS = schoolfield_original_simple(dataset, i)
        print(model_OS)
        model_OS.plot()
        if dataset.shape[0] > 5: #This model has two additional variables
            model_O = schoolfield_original(dataset, i)
            print(model_O)
            model_O.plot()