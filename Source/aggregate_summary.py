import pandas as pd
import numpy as np

print('\n\tAggregating Summary...\n')

#Geometric mean robust to 0s
def gmean(data_series, PropZero = False):
    orig_length = len(data_series)
    data_series = np.array([i for i in data_series if i != 'NA'])
    
    if PropZero:
        if not numpy.all(data_series):
            return 0
    else:
        data_series = data_series[np.nonzero(data_series)]
    
    result = np.exp(sum(np.log(data_series) / orig_length))
    return result

def mesophile_thermophile(entry):
    if entry > 323.15:
        return 'Thermophile'
    else:
        return 'Mesophile'    
    
agg_data = pd.read_csv('..\Data\summaries\summary.csv', encoding = "ISO-8859-1", index_col=0)
colnames = list(agg_data.columns.values)
#columns to keep
keep = ['Species', 'Trait', 'ConKingdom', 'ConPhylum', 'ConClass', 'ConOrder', 'ConFamily', 'ConGenus', 'Best_Guess', 'Max.response', 'Est.Tpk', 'Est.Tmin', 'Est.Tmax', 'Rank']

#Drop any column not in the keep list
for colname in colnames:
    if colname not in keep:
        agg_data.drop(colname, 1, inplace=True)
       
#select best fits for growth rates only
best = agg_data[(agg_data.Rank == 1) & (agg_data.Trait == 'Specific Growth Rate')]
#Remove blanks in best guess (metabolism)
best['Best_Guess'] = list(best['Best_Guess'].fillna('Unknown'))
print(list(best['Best_Guess']))
groups = best.groupby(['Species', 'Trait', 'ConKingdom', 'ConPhylum', 'ConClass', 'ConOrder', 'ConFamily', 'ConGenus', 'Best_Guess'])
groups = groups.aggregate(gmean)

#sort alphabetically by species name
groups.index.name = 'Species'
groups.sort_index()

#Add a temperature preference column
groups['TempPref'] =  groups['Est.Tpk'].apply(mesophile_thermophile)

#Fill blanks
groups = groups.fillna('NA')

#Save csv
groups.to_csv('../Data/summaries/aggregate_data.csv')   