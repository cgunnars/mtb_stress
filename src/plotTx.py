from util import io_feat, parse_list
import pandas as pd
import numpy as np
import seaborn as sns; sns.set()

args = io_feat()
ln_tpm = pd.read_csv(args.i, sep='\t', index_col=0)
ln_tpm.index = ln_tpm.index.str.lower()
feat_list = [x.lower() for x in parse_list(args.f)]#.lower()
print(feat_list)
#current bug is features are lowercase
ln_tpm_feat = ln_tpm.loc[feat_list, :].dropna(axis=0, how='all')

fig = sns.clustermap(data=ln_tpm_feat.T, z_score=0, xticklabels=True, yticklabels=True)
fig.savefig(args.o, bbox_inches='tight')


