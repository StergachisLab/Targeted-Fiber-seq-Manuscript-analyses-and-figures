import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load data. Firs line is the header
data_large = pd.read_csv("/mmfs1/gscratch/stergachislab/bohaczuk/scripts/publication_repositories/Targeted-Fiber-seq/DM1_CpG/largeRE_haptagged.tsv", sep="\t", header=0)
data_small = pd.read_csv("/mmfs1/gscratch/stergachislab/bohaczuk/scripts/publication_repositories/Targeted-Fiber-seq/DM1_CpG/smallRE_haptagged.tsv", sep="\t", header=0)


# Calculate -log10(p-value) and save as a new column
data_large['log_pval_large'] = -np.log10(data_large['asm_fishers_pvalue'])
data_small['log_pval_small'] = -np.log10(data_small['asm_fishers_pvalue'])




# Plot volcano plot
fig, ax = plt.subplots()
ax.scatter(data_large['mean_meth_delta'], data_large['log_pval_large'], color='#EC2027', label='Gen II/III Expanded')
ax.scatter(data_small['mean_meth_delta'], data_small['log_pval_small'], color='#006838', label='Gen I Expanded')
plt.axhline(y=-np.log10(0.05), color='black', linestyle='--')
plt.xlabel('mCpG difference')
# set x axis limits
#plt.xlim(-0.65, 0.65)
plt.ylabel('-log10(p-value)')
plt.legend()

# Save the plot
plt.savefig('volcano_plot.png')
plt.show()