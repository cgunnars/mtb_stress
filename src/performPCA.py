from sklearn import preprocessing, metrics, manifold, cluster, decomposition
import pandas as pd
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import numpy as np
from util import io_feat, parse_list, filter_genes

args = io_feat()
ln_tpm = pd.read_csv(args.i, sep='\t', index_col=0)
feat_list = parse_list(args.f)
ln_tpm = ln_tpm.loc[feat_list, :].dropna(axis=0, how='all')
z_ln_tpm = preprocessing.scale(ln_tpm.T, axis=1)
pca = decomposition.PCA()
X = pca.fit(z_ln_tpm).transform(z_ln_tpm)
label = [x[0:-2] for x in ln_tpm.columns]
# get numerical encoding for unique values
unique, c = np.unique(label, return_inverse=True)

#Plot points in PC space
ax = sns.scatterplot(x=X[:,0], y=X[:,1], hue=label, palette = 'husl')
lgd = plt.legend(bbox_to_anchor=(1.05, 1), ncol=2)
plt.xlabel("PC 1: (%.1f %% variance)" %(pca.explained_variance_ratio_[0] * 100))
plt.ylabel("PC 2: (%.1f %% variance)" %(pca.explained_variance_ratio_[1] * 100))
plt.savefig(''.join([args.o, '.pdf']), bbox_inches='tight')

#Plot loadings of first two PCs
loadings = np.transpose(pca.components_) * np.sqrt(pca.explained_variance_)
f, ax = plt.subplots(2, 1, sharex=False, sharey=True, figsize=(35,10))
for i in range(len(ax)):
    loading = loadings[:, i]
    order = np.argsort(loading)
    features = np.asarray(feat_list)[order]
    loading = loading[order]
    ax[i].bar(features, loading, color=sns.color_palette()[i], label=('PC ' + str(i + 1)))
    ax[i].set_xticklabels(features, rotation='vertical')
    ax[i].legend(loc = 2)    

plt.savefig('_'.join([args.o, 'pca-loadings.pdf']), bbox_inches='tight', pad_inches=0.5)
