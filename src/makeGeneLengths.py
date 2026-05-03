import pandas
from util import std_io

args = std_io()
df = pandas.read_csv(args.i, sep='\t')
cds_df = df[df['Feature'] == 'CDS']
cds_df['Length'] = cds_df['Stop'] - cds_df['Start'] + 1
cds_df = cds_df.sort_values(by=['Locus'])
cds_df.to_csv(args.o)
