import pandas as pd;

gff = pd.read_csv('./fastq/Mycobacterium_tuberculosis_H37Rv_gff_v3.gff', sep='\t')

#eighth column is attributes, separated by Name=; Locus=;
attributes = gff.iloc[:, 8].str.split(pat=';', expand=True).iloc[:, 1].str.split(pat='=', expand=True).iloc[:, 1]

df = pd.DataFrame()
df['GeneID'] = attributes
df['Chr'] = gff.iloc[:, 0]
df['Start'] = gff.iloc[:, 3]
df['End'] = gff.iloc[:, 4]
df['Strand'] = gff.iloc[:, 6]
print(df.head())
df.to_csv('./fastq/H37Rv_saf.gff', sep='\t', index=False)
