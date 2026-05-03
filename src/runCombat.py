from util import io_feat, filter_genes
import pandas as pd
from combat import combat
import numpy as np

args = io_feat()

pheno = pd.read_table(args.f, sep='\t', index_col=0)
dat = pd.read_table(args.i, sep='\t', index_col=0)

ebat = combat(dat, pheno['batch'])
ebat.to_csv(args.o, sep="\t")

