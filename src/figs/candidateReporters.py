import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style('white')

df = pd.read_csv('./data/all_ln-filt-tpm.txt', sep='\t', index_col=0)


genes_use = ['Rv0575c', 'Rv1375','Rv1403c', 'Rv1405c', 'Rv2390c']

df = df.loc[genes_use, :]
df = pd.melt(df.reset_index(), id_vars = 'gene', value_vars=df.columns)
df.rename(columns = {'variable': 'condition', 'value': 'ln(TPM+1)'}, inplace=True)
df['condition'] = [c[0:-2] if c[-2]=='_' else c[0:-1] for c in df['condition']]


conds_use = ['37Rv_TOL', '37Rv_LoFe_1W', '37Rv_TOL_LpH', 'Extra cellular control', 'wt_6.6', 'wt_4.5']
palette = ['#ffd582', '#ffd582','#ffaa00', '#9edfff', '#9edfff', '#00a8fc']
df = df[df['condition'].isin(conds_use)]
ax = sns.stripplot(data = df, x='gene', hue='condition', y='ln(TPM+1)', dodge=True, hue_order=conds_use, palette=palette, edgecolor='black', linewidth=1)
sns.boxplot(data = df, x='gene', hue='condition', y='ln(TPM+1)', dodge=True, showbox=False, showcaps=False, whiskerprops={'lw': 0}, hue_order=conds_use, palette=palette)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[0:len(conds_use)], conds_use, bbox_to_anchor=(1.05, 1), loc='upper left')
ax.figure.savefig('fig/acid-gene-expr.png', bbox_inches='tight')
