import pandas as pd
import numpy as np

print('\n\tBuilding final database...\n')

def to_lower(d):
    return dict((k.lower(), v) for k, v in d.items())

summary = pd.read_csv('../data/summaries/summary.csv', encoding = "ISO-8859-1", na_filter=False, index_col=0)

levels = ['ConKingdom', 'ConPhylum', 'ConClass', 'ConOrder', 'ConFamily', 'ConGenus', 'Best_Guess']
path = '../Results/Maxima_fits/{}_summary.csv'

metabolism_sep = pd.read_csv('../Results/Maxima_fits/metabolism_summary.csv', encoding = "ISO-8859-1")
metabolism_sep_gp = metabolism_sep[(metabolism_sep.Model_name == 'Boltzmann Arrhenius')]
metabolism_sep_gp.loc[metabolism_sep_gp['Species'] == 'Photosystemii', 'Species'] = 'PhotosystemII'
metabolism_sepE = metabolism_sep_gp[['Species', 'ConKingdom', 'E', 'E.max', 'E.min']].set_index(['Species', 'ConKingdom']).to_dict()

csvs = [pd.read_csv(path.format(level), encoding = "ISO-8859-1") for level in levels] 
datasets = [dataset[(dataset.Model_name == 'Boltzmann Arrhenius')] for dataset in csvs]
dicts = [dataset[['Species', 'E', 'E.max', 'E.min']].set_index(['Species']).to_dict() for dataset in datasets]

E_dict = metabolism_sepE['E']
Emin_dict = metabolism_sepE['E.min']
Emax_dict = metabolism_sepE['E.max']

summary['Metabolism_Seperate_Activation'] = [E_dict.get(i, 'NA') for i in zip(summary['Best_Guess'], summary['ConKingdom'])]
summary['Metabolism_Seperate_Activation.min'] = [Emin_dict.get(i, 'NA') for i in zip(summary['Best_Guess'], summary['ConKingdom'])]
summary['Metabolism_Seperate_Activation.max'] = [Emax_dict.get(i, 'NA') for i in zip(summary['Best_Guess'], summary['ConKingdom'])]


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