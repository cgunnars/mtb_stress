from util import std_io
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import pandas as pd 

args = std_io()
df = pd.read_csv(args.i, sep='\t', index_col=0)
cv = np.nan_to_num(np.square(stats.variation(df.values, axis = 1)))
mean = np.mean(df.values, axis = 1)
#for each gene, calculate its mean & cv2
#plot log mean vs. log cv2 
plt.scatter(np.log(mean), cv, s=1)
plt.savefig(args.o)


