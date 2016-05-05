#!/usr/bin/env python3

"""
This program provides a framework for and implementation of least squares fitting of various 
thermal response models based on experimental data

Written in Python 3.5 Anaconda Distribution
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns #All sns features can be deleted without affecting program function
import os

from lmfit import minimize, Parameters, fit_report
from scipy import stats
from datetime import datetime

starttime = datetime.now()

class estimate_parameters:
    """
    This class estimates a bunch of metabolic parameters which are then passed to the models themselves. It also extracts useful data from the database which is passed to the models. 
    """
    
    #Going to define a bunch of parameters here which can then be overwritten by passing flags to __init__
    
    k = 8.62e-5 #Boltzmann constant
    Tref = 273.15 #Reference temperature - 0C
    
    Trait_Col_Name = 'StandardisedTraitName' #trait name
    X_vals_col_name = 'ConTemp' #column name to pull x values from
    Y_vals_col_name = 'StandardisedTraitValue' #column name to pull y vals from
    
    x_val_conversion = 60 * 60 * 24
    
    species_name = ''
    
    full_strain_col_name = 'Consumer'
    genus_level_col_name = 'ConGenus'
    species_level_col_name = 'ConSpecies'
    
    species_data = True
    log_data = False #Fit the model to log transformed data
    is_celcius = True #Is the input temps in celcius
    
    def __init__(self, data, index, aux_parameters_names , flags = {}):
       
        for k, v in flags.items(): #flags will overwrite the values above allowing for flexible databases
            setattr(self, k, v)
            
        self.data = self.clean_dataset(data)
       
        self.aux_parameters_names = aux_parameters_names
        
        self.get_ancillary_info()
        self.trait = self.get_single_val(self.Trait_Col_Name)
        
        self.temps = self.get_column('K') #Temperatures
        self.responses = self.get_column('Cor_Trait_Value') #Growth rates
        
        self.set_name()
        self.get_T_pk() #Assign value for T - peak
        self.estimate_T_H()
        self.estimate_T_H_L()
        self.estimate_E_init()
        self.estimate_B0()
        
    def get_single_val(self, item):
        return self.data[item][0] #Return the first item in the column
        
    def get_column(self, item):
        return np.array(self.data[item].values) #return the whole column as an np array
        
    def clean_dataset(self, data):
        "Normalise each dataset"
        #Transform temps to kelvin
        if self.is_celcius:
            data['K'] = data[self.X_vals_col_name] + 273.15
        else:
            data['K'] = data[self.X_vals_col_name]
        
        # Convert corrected value from s^-1 to d^-1
        data['Cor_Trait_Value'] = data[self.Y_vals_col_name] * self.x_val_conversion # Convert corrected value from s^-1 to d^-1
        
        #If any trait values are negative then subtract the smallest value to normalise
        minimum_trait_value  = data['Cor_Trait_Value'].min()
        
        if minimum_trait_value < 0:
            data['Cor_Trait_Value'] -= minimum_trait_value - 10E-10 #Get rid of any 0s
        else:
            data['Cor_Trait_Value'] += 10E-10 #If there are 0s the above will not work so correct anyway
            
        if self.log_data:
            data['Cor_Trait_Value'] = np.log(data['Cor_Trait_Value'])
            data['K'] = np.log(data['K'])
            
        return data
        
    def get_ancillary_info(self):
        "Get information on each curve to include in the summary"
        self.aux_parameters_values = [self.data[aux_parameter][0] for aux_parameter in self.aux_parameters_names]

    def get_T_pk(self):
        "Find the temperature at which maximum response is observed"
        self.Tpk_row = self.data['Cor_Trait_Value'].idxmax() #Index of max response
        self.T_pk = self.data.loc[self.Tpk_row]['K'] #Temperature at which rate is maximum
        
    def estimate_T_H(self):
        "Estimate the temperature at which half the enzyme is inactivated by high temperatures"
        #slice data so only the downwards slope is included (including TPK value) 
        downslope = self.data.loc[self.Tpk_row:]
        if downslope.shape[0] > 1:
            x = downslope['K']
            
            #Linearise the slop by taking the square root
            y = downslope['Cor_Trait_Value'] ** 0.5
            
            #Calculate a regression line and find the x intercept
            slope, intercept, *vals = stats.linregress(x,y) 
            x_intercept = -intercept/slope

            #find the value of K a third of the way between T growth_max and T growth_0 (high).
            self.T_H = self.T_pk + ((x_intercept - self.T_pk) / 3) 
        else:
            #Use an arbitary value and hope for the best
            self.T_H = self.T_pk + 5

    def estimate_T_H_L(self):
        "Estimate the temperature at which half the enzyme is inactivated by low temperatures"
        #slice data so only the downwards slope is included (including TPK value)
        upslope = self.data.loc[:self.Tpk_row]
        if upslope.shape[0] > 1:
            x = upslope['K']
            #Linearise using square root
            y = upslope['Cor_Trait_Value'] ** 0.5

            slope, intercept, *vals = stats.linregress(x,y) 
            x_intercept = -intercept/slope 

            #find the value of K a quarter of the way between T growth_max and T growth_0 (low)
            self.T_H_L = self.T_pk - ((self.T_pk - x_intercept) / 4) 
        else:
            #Once again, use an arbitary value and hope for the best
            self.T_H_L = self.T_pk - 5
            
    def estimate_E_init(self):
        "Estimate energy value using the slope of the values to the peak of an arrhenius plot"
        upslope = self.data.loc[:self.Tpk_row] #slice data so only values less than TKP row are included
        if upslope.shape[0] > 1:
            temps = upslope['K']
            responses = upslope['Cor_Trait_Value']
            
            x = 1 / (self.k * temps)
            y = np.log(responses)
            
            slope, intercept, r_value, p_value, stderr = stats.linregress(x,y)
            
            self.E_init = abs(slope)
        else:
            self.E_init = 0.6 #Default value
            
    def estimate_B0(self):
        "Returns the response at the tempetature closest to Tref"
        closest_T_index = abs(self.temps - self.Tref).argmin()
        self.B0 = np.log(self.responses[closest_T_index])
            
    def set_name(self):
        "Set species name to be applied to plot title"
        if self.species_name == '' and isinstance(self.species_name, str):
            genus = self.get_single_val(self.genus_level_col_name)
            species = self.get_single_val(self.species_level_col_name)
            consumer = self.get_single_val(self.full_strain_col_name)
        
            #Use this to remove pseudoreplicates
            if pd.isnull(species) or not self.species_data:
                self.species_name = consumer #if no species is available we have to use consumer
            else:
                self.species_name = ' '.join([genus, species])
            try:    
                self.species_name = self.species_name[0].upper() + self.species_name[1:].lower() #Usually a genus and species name so this should be correct in most cases
            except TypeError:
                print('Warning, no name found at this level for group')
            
    def __str__(self):
        vars = [self.species_name, self.temps, self.responses, self.trait, self.B0, self.E_init, self.T_H_L, self.T_H, self.T_pk, self.Tpk_row]
        
        text = """\
        ----------------------
        {0[0]}
        Trait: {0[3]}
        
        Temps:
        {0[1]}
        Responses:
        {0[2]}
        
        B0: {0[4]}
        E: {0[5]}
        THL: {0[6]}
        TH: {0[7]}
        TPK: {0[8]}
        TPK Row: {0[9]}
        """.format(vars)
        
        return text
        
class physiological_growth_model:
    k = 8.62e-5 #Boltzmann constant
    Tref = 273.15 #Reference temperature - 0C
    
    response_corrected = False #In cases where a peak in the curve tends to infinity this will be set to true. 
    rank = 'NA' #When fitting multiple models to the same data we will use this to rank models
    
    def extract_parameters(self, est_parameters):
        self.temps = est_parameters.temps
        self.responses = est_parameters.responses
        self.species_name = est_parameters.species_name
        self.trait = est_parameters.trait
        self.B0 = est_parameters.B0
        self.E_init = est_parameters.E_init
        self.T_H_L = est_parameters.T_H_L
        self.T_H = est_parameters.T_H
        self.T_pk = est_parameters.T_pk
        self.Tpk_row = est_parameters.Tpk_row
        self.aux_parameters_names = est_parameters.aux_parameters_names
        self.aux_parameters_values = est_parameters.aux_parameters_values
         
    def get_residuals(self, params, temps, responses):
        "Called by fit model only, generates residuals using test values"
        residuals = np.exp(self.fit(params, self.temps)) - responses
        
        return residuals
        
    def fit_model(self):
        "Least squares regression to minimise fit"
        self.model = minimize(self.get_residuals, 
                              self.parameters, 
                              args=(self.temps, self.responses),
                              method="leastsq")
                              
        self.R2 = 1 - np.var(self.model.residual) / np.var(self.responses)
        
    def assess_model(self):
        """Calculate the Akaike Information Criterion and Bayesian Information Criterion, using:

        - n:   number of observations
        - k:   number of parameters
        - rss: residual sum of squares
        """
        k = self.model.nvarys #Number of variables
        n = self.model.ndata #Number of data points
        rss = sum(np.power(self.model.residual, 2)) #Residual sum of squares
        
        self.AIC = 2 * k + n * np.log(rss / n)
        self.BIC = np.log(n)*(k + 1) + n * np.log(rss / n)
           
    def smooth(self):
        "Pass an interpolated list of temperature values back through the curve function to generate a smooth curve"
        self.smooth_x = np.arange(self.temps.min() - 3, self.temps.max() + 3, 0.1) #Extrapolate a little 
        self.smooth_y = np.exp(self.fit(self.model.params, self.smooth_x))
        
    def est_final_temps(self):
        """
        Slightly messy function as schoolfield cannot be evalutated non numerically
        First finds an estimate of TPK
        Uses the TPK value to split the growth curve into an upwards slope and a downwards slope
        Then finds the temperature at the rate which is the given percentage of the maximum rate
        """ 
        percentile = 0.75 #Closest percentile to find 
        
        #check 15 degrees either side, if its not in this range it won't be found anyway
        peak_check_x = np.arange(self.temps.min() - 15, self.temps.max() + 15, 0.1)
        peak_check_y = np.exp(self.fit(self.model.params, peak_check_x))

        #Check the max response isn't the final value in the curve
        if peak_check_y.argmax() < peak_check_y.size - 1:
            max_index = peak_check_y.argmax()
            max_val = peak_check_y.max()
            
            #Index the temperatures with the maximum response value's position
            self.tpk_est = peak_check_x[max_index]    
            self.max_response_est = max_val
            
            #In some cases the fitting function leads to a near infinitisimal peak being found skewing the value
            #Check if that has happened and correct to the maximum observed rate if it has!
            
            max_val_check_1 = peak_check_y[max_index - 5] / self.max_response_est
            
            #Kicks in if response increases 15% in half a degrees
            if max_val_check_1 < 0.85:
                self.response_corrected = True
                self.max_response_est = self.responses.max()

            #Find each y value as a percentage of the maximum
            upslope_y = peak_check_y[:max_index + 1] / max_val
            downslope_y = peak_check_y[max_index:] / max_val
            upslope_x = peak_check_x[:max_index + 1]
            downslope_x = peak_check_x[max_index:]
            
            #Get the closest value to the percentage
            position_lower = (np.abs(upslope_y-percentile)).argmin()
            position_upper = (np.abs(downslope_y-percentile)).argmin()
            
            #Index the original x values to find percentile
            self.lower_percentile = upslope_x[position_lower]
            self.upper_percentile = downslope_x[position_upper]
            
        else:
            #Set default values
            self.tpk_est = 'NA'
            self.lower_percentile = 'NA'
            self.upper_percentile = 'NA'
            self.max_response_est = 'NA'
            
    def plot(self, out_path, scale_type='standard'):
        #apologies to anyonne reading this function
        #We can use this to ID which model is best
        textdata = [self.final_E, self.final_B0, self.R2, self.AIC, self.BIC] #Added to plot to show fit quality
        title = '{}: {}'.format(self.index, self.species_name) 
        
        plt_x = self.temps
        plt_y = self.responses
        
        plt_x_curve = self.smooth_x
        plt_y_curve = self.smooth_y
        
        if scale_type == 'log':
            plt_x, plt_y, plt_x_curve, plt_y_curve = np.log(plt_x), np.log(plt_y), np.log(plt_x_curve), np.log(plt_y_curve)
                
        if scale_type == 'arrhenius':
            plt_y, plt_y_curve = np.log(plt_y), np.log(plt_y_curve)
            plt_x, plt_x_curve = 1 / plt_x, 1 / plt_x_curve
        
        f = plt.figure()
        sns.set_style("ticks", {'axes.grid': True})
        ax = f.add_subplot(111)
        
        #Plot actual observed data as a scatterplot
        plt.plot(plt_x, plt_y, marker='o', linestyle='None', color='green', alpha=0.7)
        
        #Plot fitted curve
        plt.plot(plt_x_curve, plt_y_curve, marker='None', color='royalblue', linewidth=3)
        
        plt.title(title, fontsize=14, fontweight='bold')
        
        if scale_type == 'log':
            plt.xlabel('log(Temperature (K))')
            plt.ylabel('Log(' + self.trait + ')')
                
        if scale_type == 'arrhenius':
            plt.xlabel('1 / Temperature (K)')
            plt.ylabel('Log(' + self.trait + ')')
            
        else:
            plt.xlabel('Temperature (K)')
            plt.ylabel(self.trait)
        
        plt.text(0.05, 0.85,'E:  {0[0]:.2f}\nB0:  {0[1]:.2f}\nR2:  {0[2]:.2f}\nAIC: {0[3]:.2f}\nBIC: {0[4]:.2f}'.format(textdata),
                 ha='left', va='center', transform=ax.transAxes, color='darkslategrey')
        sns.despine() #Remove top and right border
        
        output_path = out_path + '/{}.png'.format(self.index)
        plt.savefig(output_path, bbox_inches='tight') 
        plt.close()
        
    def get_stderrs(self):
        "These aren't actually outputted anywhere, but it would be easy enough to make them so I'm leaving this here"
        self.final_B0_stderr = self.model.params['B0_start'].stderr
        self.final_E_stderr = self.model.params['E'].stderr 
            
    def __lt__(self, other):
        "By defining less than we can directly compare models using max() to determine which is best"
        if self.AIC == other.AIC:
            return self.BIC > other.BIC
        else:
            return self.AIC > other.AIC
            
    def __eq__(self, other):
        return self.AIC == other.AIC and self.BIC == other.BIC

    
class Boltzmann_Arrhenius(physiological_growth_model):
    "Simplest model - only works with linear upwards responses"
    model_name = "Boltzmann Arrhenius"
    def __init__(self, est_parameters, index):
        self.extract_parameters(est_parameters)
        self.index = str(index) + "_BA" #ID for the model
        self.set_parameters()
        self.fit_model()
        self.smooth()
        self.assess_model()
        self.get_final_values()
        self.get_stderrs()
        
    def set_parameters(self):
        "Create a parameters object using out guesses, these will then be fitted using least squares regression"
        self.parameters = Parameters()
        #                             Name,  Start,   Can_Vary, Lower, Upper
        self.parameters.add_many(('B0_start', self.B0, True, -np.inf, np.inf,  None),
                                 ('E', self.E_init, True, 0, np.inf,  None))

    def fit(self, params, temps):
        "Fit a schoolfield curve to a list of temperature values"
        parameter_vals = params.valuesdict()
        B0 = parameter_vals['B0_start'] #Basic metabolic rate
        E = parameter_vals['E'] #Activation energy of enzymes

        fit = B0 - E / self.k * (1 / temps - 1 / self.Tref)

        return fit
        
    def get_final_values(self):
        "Get the final fitted values for the model"
        values = self.model.params.valuesdict()
        self.final_B0 = values['B0_start']
        self.final_E = values['E']
        
    def __str__(self):
        "Allows print() to be called on the object"
        vars = [self.species_name, self.B0, self.final_B0, self.E_init, self.final_E, self.R2, self.AIC, self.BIC]
        text = """\
        ---Boltzmann Arrhenius Model---
        {0[0]}
        
        B0 est = {0[1]:.2f}
        B0 final = {0[2]:.2f}
        
        E est = {0[3]:.2f}
        E final = {0[4]:.2f}
        
        R2: = {0[5]:.2f}
        AIC = {0[6]:.2f}
        BIC = {0[7]:.2f}
        
        """.format(vars)
        return text        
        
class schoolfield_two_factor(physiological_growth_model):
    "Schoolfield model using T_pk as a substitute for T_H"
    model_name = "schoolfield two factor"
    def __init__(self, est_parameters, index):
        self.extract_parameters(est_parameters)
        self.index = str(index) + "_Sch_TF" #ID for the model
        self.set_parameters()
        self.fit_model()
        self.smooth()
        self.assess_model()
        self.get_final_values()
        self.get_stderrs()
        self.est_final_temps()
        
    def set_parameters(self):
        "Create a parameters object using out guesses, these will then be fitted using least squares regression"
        self.parameters = Parameters()
        #                   Name,      Start,   Can_Vary, Lower, Upper
        self.parameters.add_many(('B0_start', self.B0, True, -np.inf, np.inf,  None),
                          ('E', self.E_init, True, 10E-10, np.inf,  None),
                          ('E_D',self.E_init * 4, True, 10E-10, np.inf,  None),
                          ('T_pk', self.T_pk, True, 273.15-50, 273.15+150,  None))

    def fit(self, params, temps):
        "Fit a schoolfield curve to a list of temperature values"
        parameter_vals = params.valuesdict()
        B0 = parameter_vals['B0_start'] #Basic metabolic rate
        E = parameter_vals['E'] #Activation energy of enzymes
        E_D = parameter_vals['E_D'] #Inactivation energy of enzymes
        T_pk = parameter_vals['T_pk'] #Temperature at which peak response is observed  
        
        fit = B0 + np.log(np.exp((-E / self.k) * ((1 / temps) - (1 / self.Tref))) /\
                         (1 + (E/(E_D - E)) * np.exp(E_D / self.k * (1 / T_pk - 1 / temps))))

        return fit
        
    def get_final_values(self):
        "Get the final fitted values for the model"
        values = self.model.params.valuesdict()
        self.final_B0 = values['B0_start']
        self.final_E = values['E']
        self.final_E_D = values['E_D']
        self.final_T_pk = values['T_pk']   

    def get_stderrs(self):
        "These aren't actually outputted anywhere, but it would be easy enough to make them so I'm leaving this here"
        self.final_B0_stderr = self.model.params['B0_start'].stderr
        self.final_E_stderr = self.model.params['E'].stderr 
        self.final_E_D_stderr = self.model.params['E_D'].stderr
        self.final_T_pk_stderr = self.model.params['T_pk'].stderr
        
    def __str__(self):
        "Allows print() to be called on the object"
        vars = [self.species_name, self.B0, self.final_B0, self.E_init, self.final_E, self.T_pk, self.final_T_pk,
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
        
        """.format(vars)
        return text

class schoolfield_original_simple(physiological_growth_model):
     # Original Sharpe-Schoolfield's model with low-temp inactivation term removed, very similar to previous model
     # This seems to generate pretty much identical results to the two factor model, but throws its toys out the pram sometimes.
     # In this version T_H is the temperature at which half the enzyme is denatured by heat stress
    model_name = "schoolfield simple"
    def __init__(self, est_parameters, index):
        self.extract_parameters(est_parameters)
        self.index = str(index) + "_Sch_OS" #Used to name plot graphics file
        self.set_parameters()
        self.fit_model()
        self.smooth()
        self.assess_model()
        self.get_final_values()
        self.get_stderrs()
        self.est_final_temps()
        
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
        B0 = parameter_vals['B0_start'] #Basic metabolic rate
        E = parameter_vals['E'] #Activation energy of enzymes
        E_D = parameter_vals['E_D'] #Inactivation energy of enzymes
        T_H = parameter_vals['T_H'] #Temperature at which half od enzzymes are denatured
        
        fit = B0 + np.log(np.exp((-E / self.k) * ((1 / temps) - (1 / self.Tref)))\
                        /(1 + np.exp((E_D / self.k) * (1 / T_H - 1 / temps))))
        
        return fit
        
    def get_final_values(self):
        "Get the final fitted values for the model"
        values = self.model.params.valuesdict()
        self.final_B0 = values['B0_start']
        self.final_E = values['E']
        self.final_E_D = values['E_D']
        self.final_T_H = values['T_H']   
        
    def get_stderrs(self):
        "These aren't actually outputted anywhere, but it would be easy enough to make them so I'm leaving this here"
        self.final_B0_stderr = self.model.params['B0_start'].stderr
        self.final_E_stderr = self.model.params['E'].stderr 
        self.final_E_D_stderr = self.model.params['E_D'].stderr
        self.final_T_stderr = self.model.params['T_H'].stderr
        
    def __str__(self):
        "Allows print() to be called on the object"
        vars = [self.species_name, self.B0, self.final_B0, self.E_init, self.final_E, self.T_pk, self.T_H, 
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

class schoolfield_original(physiological_growth_model):
     # Original Sharpe-Schoolfield's model with low-temp inactivation term removed, very similar to previous model
     # In this version T_H_L is a low temperature enzyme inactivation constant (as if having a high temp one wasn't fun enough already)
    model_name = "schoolfield"
    def __init__(self, est_parameters, index):
        self.extract_parameters(est_parameters)
        self.index = str(index) + "_Sch_O" #Used to name plot graphics file
        self.set_parameters()
        self.fit_model()
        self.smooth()
        self.assess_model()
        self.get_final_values()
        self.get_stderrs()
        self.est_final_temps()
        
    def set_parameters(self):
        "Create a parameters object using out guesses, these will then be fitted using least squares regression, note additional T_H_L parameter"
        self.parameters = Parameters()
        #                   Name,      Start,   Can_Vary, Lower, Upper
        self.parameters.add_many(('B0_start', self.B0, True, -np.inf, np.inf,  None),
                          ('E', self.E_init, True, 10E-10, np.inf,  None),
                          ('E_D',self.E_init * 4, True, 10E-10, np.inf,  None),
                          ('E_D_L',self.E_init * (-6), True, -np.inf, -10E-10,  None),
                          ('T_H', self.T_H, True, self.T_pk + 0.1, 273.15+170,  None),
                          ('T_H_L', self.T_H_L, True, 273.15-70, self.T_pk - 0.1,  None))

    def fit(self, params, temps):
        "Fit a schoolfield curve to a list of temperature values"
        parameter_vals = params.valuesdict()
        B0 = parameter_vals['B0_start'] #Basic metabolic rate
        E = parameter_vals['E'] #Activation energy of enzymes
        E_D = parameter_vals['E_D'] #Inactivation energy of enzymes
        E_D_L = parameter_vals['E_D_L'] #Energy of cold inactivation
        T_H = parameter_vals['T_H'] #Temperature at which half of enzymes are denatured
        T_H_L = parameter_vals['T_H_L'] #Temperature at which half of enzymes are cold inactivated
        
        fit = B0 + np.log(np.exp((-E / self.k) * ((1 / temps) - (1 / self.Tref)))    \
                        /(1 + np.exp((E_D_L / self.k) * (1 / T_H_L - 1 / temps)) +   \
                              np.exp((E_D / self.k) * (1 / T_H - 1 / temps))))
        
        return fit
        
    def get_final_values(self):
        "Get the final fitted values for the model"
        values = self.model.params.valuesdict()
        self.final_B0 = values['B0_start']
        self.final_E = values['E']
        self.final_E_D = values['E_D']
        self.final_E_D_L = values['E_D_L']
        self.final_T_H = values['T_H']
        self.final_T_H_L = values['T_H_L']

    def get_stderrs(self):
        "These aren't actually outputted anywhere, but it would be easy enough to make them so I'm leaving this here"
        self.final_B0_stderr = self.model.params['B0_start'].stderr
        self.final_E_stderr = self.model.params['E'].stderr 
        self.final_E_D_stderr = self.model.params['E_D'].stderr
        self.final_E_D_L_stderr = self.model.params['E_D_L'].stderr
        self.final_T_stderr = self.model.params['T_H'].stderr        
        self.final_T_H_L_stderr = self.model.params['T_H_L'].stderr    
        
    def __str__(self):
        "Allows print() to be called on the object"
        vars = [self.species_name, self.B0, self.final_B0, self.E_init, self.final_E, self.E_init * 4, self.final_E_D, self.E_init * (-6), self.final_E_D_L,
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
            
def get_datasets(path, sep = 'OriginalID', _sort = ['ConTemp']):
    "Create a set of temperature response curve datasets from csv"
    data = pd.read_csv(path, encoding = "ISO-8859-1") #Open in latin 1
    
    data['FinalID'] = pd.factorize(data[sep])[0]
    ids = pd.unique(data['FinalID']).tolist() #Get unique identifiers
    
    #create a dictionary of datasets for easy access later
    Datasets = {}
    for id in ids:
        curve_data = data.loc[data['FinalID'] == id] #seperate data by uniqueID
        curve_data = curve_data.sort_values(_sort).reset_index() #sort so rows are in temperature order, reset index to 0  
        Datasets[id] = curve_data
    return Datasets    

def rank_and_flatten(model_list):
    "Function to rank a nested lit of models and flatten the list"
    all_models = []
    if any(isinstance(i, list) for i in model_list): #check if list is nested
        for candidates in model_list:
            candidates.sort() #Uses built in lt method, best model will be first
            for rank, model in enumerate(candidates[::-1]):
                model.rank = rank + 1 #setting post hoc
                all_models.append(model)
    else:
        all_models = model_list
    return all_models
    
def output_csv(model_list, aux_cols = [], path = 'summary.csv', whole_curves = False, sortby=['Species', 'FinalID']):
    col_names = ["Species", "Model_name", "Trait", "B0", "E", "T_pk", "E_D", "E_D_L", "Est.Tpk",
                 "Est.Tmin", "Est.Tmax", "Max response", "R_Squared", "AIC", "BIC",
                 "Rank", "Corrected", "Number.of.Data.Points", "Number.of.Variables"] + aux_cols
                 
    if whole_curves:
        col_names += ['Temperature', 'Response', 'Original_Data']
    
    rows = []
    
    for model in model_list:
        final_T_pk = getattr(model, 'final_T_pk', "NA")
        final_estimated_T_pk = getattr(model, 'tpk_est', "NA")
        final_max_response = getattr(model, 'max_response_est', "NA")
        final_upper_percentile = getattr(model, 'upper_percentile', "NA")
        final_lower_percentile = getattr(model, 'lower_percentile', "NA")
        final_E_D = getattr(model, 'final_E_D', "NA")
        final_E_D_L = getattr(model, 'final_E_D_L', "NA")
        final_T_H = getattr(model, 'final_T_H', "NA")
        final_T_H_L = getattr(model, 'final_T_H_L', "NA")
        
        if whole_curves: #output the entire smooth growth curve
            for i, x in np.ndenumerate(model.smooth_x):
                temp = x
                resp = model.smooth_y[i]
                key = model.index + '_model_point_' + str(i[0]) #needs to be unique
                
                entries = [model.species_name, model.model_name, model.trait, model.final_B0, model.final_E, final_T_pk, final_E_D, final_E_D_L, final_estimated_T_pk, 
                final_lower_percentile, final_upper_percentile, final_max_response, temp, resp, model.R2, model.AIC, model.BIC, model.rank,
                model.response_corrected, model.model.ndata, model.model.nvarys] + model.aux_parameters_values + [temp, resp, 'False']
                
                row = (key, entries)
                rows.append(row)

            for i, x in np.ndenumerate(model.temps):
                temp = x
                resp = model.responses[i]
                key = model.index + '_orig_point_' + str(i[0]) #needs to be unique
                
                entries = [model.species_name, model.model_name, model.trait, model.final_B0, model.final_E, final_T_pk, final_E_D, final_E_D_L, final_estimated_T_pk, 
                final_lower_percentile, final_upper_percentile, final_max_response, temp, resp, model.R2, model.AIC, model.BIC, model.rank,
                model.response_corrected, model.model.ndata, model.model.nvarys] + model.aux_parameters_values + [temp, resp, 'True']
                
                row = (key, entries)
                rows.append(row)   
                
        else:
            entries = [model.species_name, model.model_name, model.trait, model.final_B0, model.final_E, final_T_pk, final_E_D, final_E_D_L, 
                       final_estimated_T_pk, final_lower_percentile, final_upper_percentile, final_max_response, model.R2, model.AIC, model.BIC, 
                       model.rank, model.response_corrected, model.model.ndata, model.model.nvarys] + model.aux_parameters_values 
                   
            row = (model.index, entries)
            rows.append(row)
        
    df = pd.DataFrame.from_items(rows,  orient='index', columns=col_names)
    df = df.sort_values(sortby)
    df.to_csv(path)    