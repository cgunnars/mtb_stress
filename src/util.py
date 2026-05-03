def std_io():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required = True)
    parser.add_argument('-o', required = True)
    args = parser.parse_args()
    return args

#TODO base fxn + wrappers
def io_feat():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required = True)
    parser.add_argument('-f', required = True)
    parser.add_argument('-o', required = True)
    return parser.parse_args()

def parse_list(f):
    with open(f) as f:
        lines = [line.strip() for line in f]
    return lines

def filter_genes(df, thresh=25):
    import numpy as np
    import pandas as pd
    print(df.head())
    keep = np.sum(df.values >= thresh, axis = 1) >= 2
    df = df.loc[keep, :]
    return df

#return top N most variable genes by CV2
def get_var_genes(df, thresh = 500):
    import scipy.stats as stats
    import numpy as np
    # need to deal with mean = 0, std = 0 cases
    cv = np.nan_to_num(np.square(stats.variation(df.values, axis = 1)))
    order = np.argsort(cv)
    topn = order[-(thresh + 1):-1]
    return df.index[topn]


