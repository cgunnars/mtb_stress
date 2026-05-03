import argparse
import os
import pandas as pd
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs='+')
parser.add_argument('-o')

args = parser.parse_args()

genes_df = pd.DataFrame()
genes = []
for filename in args.i:
    splits = os.path.basename(filename).split('_')
    exp = splits[0]
    cond = splits[1]
    spec = splits[2]
    with open(filename, 'r') as f:
        genes_f = [line.strip().lower() for line in f] 
        genes = list(set(genes + genes_f))
        f_df = pd.DataFrame(data = np.repeat(True, len(genes_f)), index=genes_f, columns=[filename])
        genes_df = pd.concat([genes_df, f_df], axis=1, sort=True) 
    

genes_df = genes_df.fillna(False)
genes_df['Total'] = genes_df.sum(axis=1)
genes_df.to_csv(args.o + '.csv')
print(genes_df.shape)
if (genes_df.shape[1] <= 3):
    print('small')
    genes_subset = genes_df[genes_df['Total'] >= 2]
# expressed in all but two conditions
else:
    genes_subset = genes_df[genes_df['Total'] >= (genes_df.shape[1] - 3)]
with open(args.o + '.txt', 'w') as f:
    f.writelines('\n'.join(genes_subset.index))


