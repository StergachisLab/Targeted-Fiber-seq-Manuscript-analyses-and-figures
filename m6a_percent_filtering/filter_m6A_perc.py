import pandas as pd
import numpy as np
import sys

filename=sys.argv[1]
low_bound=sys.argv[2]
up_bound=sys.argv[3]

data=pd.read_csv(filename,sep='\t')

data['m6Afrac']=data['total_m6a_bp']/data['total_AT_bp']


# Make a list of reads to include
include_reads = data.loc[(data['m6Afrac'] >= float(low_bound)) & (data['m6Afrac'] <= float(up_bound)), 'read_name'].tolist()
print('include',include_reads)

# write list to file
with open('include_list.txt','w') as f:
    for read in include_reads:
            f.write(read+'\n')

