import pandas 
import argparse


def countToTpm(mat, feat, mean_rl = 100, locus=False):
    if locus:
        txlookup = 'Locus'
    else:
        txlookup = 'Name'
    validtx = [tx for tx in feat[txlookup] if tx in mat.index.tolist()]
    feat = feat[feat[txlookup].isin(validtx)]
    mat = mat[mat.index.isin(validtx)]
    
    rpk = mat.div(feat['Length'].values, axis = 0)
    total_count = rpk.sum(axis = 0, numeric_only=True) / (10 ** 6)
    tpm = rpk.div(total_count.values.T, axis = 1) #* mean_rl * (10 ** 6)
    if locus: #set index to gene name so consistent
        tpm.index = feat['Name']
        

    return tpm
    
parser = argparse.ArgumentParser()
parser.add_argument('-i', help='counts file')
parser.add_argument('-g', help='gene lengths file')
parser.add_argument('--locus', action='store_true')
parser.add_argument('-o', help='output')
args = parser.parse_args()

counts = pandas.read_csv(args.i, sep='\t', index_col=0)
print(counts.head())
features = pandas.read_csv(args.g)
tpm = countToTpm(counts, features, locus=args.locus)
tpm.to_csv(args.o, sep='\t')
