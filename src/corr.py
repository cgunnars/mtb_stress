from util import io_feat, parse_list
from sklearn import preprocessing, metrics, manifold, cluster, decomposition
import pandas as pd
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import numpy as np

args = io_feat()
ln_tpm = pd.read_csv(args.i, sep='\t', index_col=0)
feat_list = parse_list(args.f) 
ln_tpm = ln_tpm.loc[feat_list,:]
z_ln_tpm = pd.DataFrame(data=preprocessing.scale(ln_tpm, axis=0), index=ln_tpm.index, columns=ln_tpm.columns)
corr = z_ln_tpm.corr(method='pearson')
fig = plt.figure(figsize=(15,13))
fig = sns.clustermap(data=corr, square=True, cmap='Blues', cbar_kws={"shrink": .5},
                 xticklabels=True, yticklabels=True)
fig.savefig(args.o, bbox_inches='tight', figsize = (40,40))

