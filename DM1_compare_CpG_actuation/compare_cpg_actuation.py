import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy

# Load CpG data. Firs line is the header
data_cpg = pd.read_csv("largeRE_haptagged.tsv", sep="\t", header=0)

# rename chrom to chr
data_cpg = data_cpg.rename(columns={'chrom':'chr'})


# Exclude the first row which straddles a cut site, and reindex the dataframe
data_cpg = data_cpg.iloc[1:]
data_cpg = data_cpg.reset_index(drop=True)



# Load long RE actuation data. First line is the header
data_actuation = pd.read_csv("../compare_perc_actuation/DM1/peaks/normal_peaks_widepeak_stats_with_pvals.txt", sep="\t", header=0)

# merge the two dataframes
data = pd.merge(data_cpg, data_actuation, on=['chr', 'start', 'end'], how='inner')


data.loc[data["start"] == 45769240, 'label_large'] = 'SIX5 promotor'
data.loc[data["start"] == 45770259, 'label_large'] = "CTCF peak"
data.loc[data["start"] == 45770597, 'label_large'] = "putative SIX5 enhancer"
data.loc[data["start"] == 45810188, 'label_large'] = 'RSPA6A intron 2'



# plot cpg vs actuation
fig, ax = plt.subplots()
ax.scatter(data['fraction_actuation_diff_large'], data['mean_meth_delta'], color='#EAA4AA', label='Gen II/III Expanded')
# Color by statistical significance and add labels
for i in range(len(data)):
	if data['corrected_pvals_large'][i] < 0.05:
		ax.scatter(data['fraction_actuation_diff_large'][i], data['mean_meth_delta'][i], color='#EC2027')
	if not pd.isnull(data['label_large'][i]):
		ax.text(data['fraction_actuation_diff_large'][i], data['mean_meth_delta'][i], data['label_large'][i], fontsize=8)
plt.axhline(y=0, color='black', linestyle='--')
plt.axvline(x=0, color='black', linestyle='--')
plt.ylim(-0.5,0.5)
plt.xlim(-0.65,0.65)

plt.ylabel('mCpG difference')
plt.xlabel('Fraction actuation difference')


plt.savefig('actuation_vs_cpg_large.png')
plt.savefig('actuation_vs_cpg_large.pdf', format='pdf', bbox_inches='tight')
plt.show()