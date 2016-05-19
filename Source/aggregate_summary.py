import pandas as pd
import numpy as np

print('\n\tAggregating Summary...\n')

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

agg_data = pd.read_csv('..\Data\summaries\summary.csv', encoding = "ISO-8859-1", index_col=0)
colnames = list(agg_data.columns.values)
keep = ['Species', 'Trait', 'ConKingdom', 'ConPhylum', 'ConClass', 'ConOrder', 'ConFamily', 'ConGenus', 'Best_Guess', 'Max.response', 'Est.Tpk', 'Est.Tmin', 'Est.Tmax', 'Rank']

for colname in colnames:
    if colname not in keep:
        agg_data.drop(colname, 1, inplace=True)
       
best = agg_data[(agg_data.Rank == 1) & (agg_data.Trait == 'Specific Growth Rate')]
groups = best.groupby(['Species', 'Trait', 'ConKingdom', 'ConPhylum', 'ConClass', 'ConOrder', 'ConFamily', 'ConGenus', 'Best_Guess'])
groups = groups.aggregate(gmean)
groups.index.name = 'Species'
groups.sort_index()
groups = groups.fillna('NA')

groups.to_csv('../Data/summaries/aggregate_data.csv')   