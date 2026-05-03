from util import std_io, filter_genes
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import pandas as pd
import numpy as np

args = std_io()
fig, ax = plt.subplots(1, 3, figsize=(10,3)) 

tpm = pd.read_csv(args.i, sep='\t', index_col=0)
filt_tpm = filter_genes(tpm).values.flatten()
tpm = tpm.values.flatten()
sns.distplot(tpm, kde=False, ax = ax[0], axlabel = 'TPM')
sns.distplot(np.log(tpm + 1), kde=False, ax = ax[1], axlabel = 'ln(TPM + 1)')
sns.distplot(np.log(filt_tpm + 1), kde=False, ax=ax[2], axlabel = 'ln(TPM + 1), filtered')
fig.savefig(args.o, bbox_inches='tight')

