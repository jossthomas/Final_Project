import pandas as pd
import numpy as np

print('\n\tBuilding final database...\n')

#make all items in a dictionary lower case
def to_lower(d):
    return dict((k.lower(), v) for k, v in d.items())

def mesophile_thermophile(entry):
    try:
        entry = float(entry)
        if entry > 323.15:
            return 'Thermophile'
        else:
            return 'Mesophile'
    except:
        return 'NA'       

summary = pd.read_csv('../data/summaries/summary.csv', encoding = "ISO-8859-1", na_filter=False, index_col=0)

#create a temperature preference column
summary['TempPref'] =  summary['Est.Tpk'].apply(mesophile_thermophile)

levels = ['ConKingdom', 'ConPhylum', 'ConClass', 'ConOrder', 'ConFamily', 'ConGenus']
#output path, will be string formatted later (should be updated to use a path object in future)
path = '../Results/Maxima_fits/{}_summary.csv'

#Metabolism level activations for bacteria and archaea

metabolism_sep = pd.read_csv('../Results/Maxima_fits/metabolism_summary.csv', encoding = "ISO-8859-1")
#Remove the linear models
metabolism_sep_gp = metabolism_sep[(metabolism_sep.Model_name == 'Boltzmann Arrhenius')]
#Make sure string formatting matches up
metabolism_sep_gp.loc[metabolism_sep_gp['Species'] == 'Photosystemii', 'Species'] = 'PhotosystemII'
#Create a nested dictionary containing E, E.max and E.min
metabolism_sepE = metabolism_sep_gp[['Species', 'ConKingdom', 'E', 'E.max', 'E.min']].set_index(['Species', 'ConKingdom']).to_dict()

#Extract nested dictionaries
E_dict_met = metabolism_sepE['E']
Emin_dict_met = metabolism_sepE['E.min']
Emax_dict_met = metabolism_sepE['E.max']

#Perform what is essentially a vlookup, use zip so we can chec both temp pref and kingdom
summary['Metabolism_Seperate_Activation'] = [E_dict_met.get(i, 'NA') for i in zip(summary['Best_Guess'], summary['ConKingdom'])]
summary['Metabolism_Seperate_Activation.min'] = [Emin_dict_met.get(i, 'NA') for i in zip(summary['Best_Guess'], summary['ConKingdom'])]
summary['Metabolism_Seperate_Activation.max'] = [Emax_dict_met.get(i, 'NA') for i in zip(summary['Best_Guess'], summary['ConKingdom'])]

#Mesophile thermophile level activations for bacteria and archaea, same as above

Temp_Groups = pd.read_csv('../Results/Maxima_fits/Temp_Group_summary.csv', encoding = "ISO-8859-1")
Temp_Groups = Temp_Groups[(Temp_Groups.Model_name == 'Boltzmann Arrhenius')]
Temp_dicts = Temp_Groups[['Species', 'ConKingdom', 'E', 'E.max', 'E.min']].set_index(['Species', 'ConKingdom']).to_dict()

E_dict_temp = Temp_dicts['E']
Emin_dict_temp = Temp_dicts['E.min']
Emax_dict_temp = Temp_dicts['E.max']

summary['TempPref_Activation'] = [E_dict_temp.get(i, 'NA') for i in zip(summary['TempPref'], summary['ConKingdom'])]
summary['TempPref_Activation.min'] = [Emin_dict_temp.get(i, 'NA') for i in zip(summary['TempPref'], summary['ConKingdom'])]
summary['TempPref_Activation.max'] = [Emax_dict_temp.get(i, 'NA') for i in zip(summary['TempPref'], summary['ConKingdom'])]

#Every other groups activations for bacteria and archaea, same as above but in a loop so we can do all 6 levels in one go. 

csvs = [pd.read_csv(path.format(level), encoding = "ISO-8859-1") for level in levels] 
datasets = [dataset[(dataset.Model_name == 'Boltzmann Arrhenius')] for dataset in csvs]
dicts = [dataset[['Species', 'E', 'E.max', 'E.min']].set_index(['Species']).to_dict() for dataset in datasets]

for i, level in enumerate(levels):
    colname = '{}_Activation'.format(level)
    colname_Emin = '{}.E.min'.format(level)
    colname_Emax = '{}.E.max'.format(level)
    E_dict = to_lower(dicts[i]['E'])
    Emin_dict = to_lower(dicts[i]['E.min'])
    Emax_dict = to_lower(dicts[i]['E.max'])
    summary[colname] = [E_dict.get(i.lower(), 'NA') for i in summary[level]]
    summary[colname_Emin] = [Emin_dict.get(i.lower(), 'NA') for i in summary[level]]
    summary[colname_Emax] = [Emax_dict.get(i.lower(), 'NA') for i in summary[level]]

summary.to_csv('../Data/summaries/summary_activation.csv')

print('\n\tFinished!\n')