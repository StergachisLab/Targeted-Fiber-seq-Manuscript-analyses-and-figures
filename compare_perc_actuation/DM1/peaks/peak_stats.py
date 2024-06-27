import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy



# Load data. Firs line is the header
data = pd.read_csv("normal_peaks_counts.txt", sep="\t", header=0)

# Exclude the first row which straddles a cut site, and reindex the dataframe
data = data.iloc[1:]
data = data.reset_index(drop=True)


# Calculate fraction actuation for each line and save as a new column
data['fraction_actuation_normal']= data['norm_FIRE']/(data['norm_FIRE']+data['norm_other'])
data['fraction_actuation_large'] = data['large_FIRE']/(data['large_FIRE']+data['large_other'])
data['fraction_actuation_small'] = data['small_FIRE']/(data['small_FIRE']+data['small_other'])


# Calculate the difference between the fraction actuation, save as a new column
data['fraction_actuation_diff_large'] = -data['fraction_actuation_normal'] + data['fraction_actuation_large']
data['fraction_actuation_diff_small'] = -data['fraction_actuation_normal'] + data['fraction_actuation_small']

# Calculate the p-value with Fisher's exact test  and save as a new column
data['pval_large'] = data.apply(lambda row: scipy.stats.fisher_exact([[row['norm_FIRE'], row['norm_other']], [row['large_FIRE'], row['large_other']]])[1], axis=1)
data['pval_small'] = data.apply(lambda row: scipy.stats.fisher_exact([[row['norm_FIRE'], row['norm_other']], [row['small_FIRE'], row['small_other']]])[1], axis=1)

# Calculate -log10(p-value) and save as a new column
data['log_pval_large'] = -np.log10(data['pval_large'])
data['log_pval_small'] = -np.log10(data['pval_small'])

# apply benjamini-hochberg correction
corrected_pvals_large=scipy.stats.false_discovery_control(data['pval_large'])
corrected_pvals_small=scipy.stats.false_discovery_control(data['pval_small'])

# add the corrected p-values to the dataframe
data['corrected_pvals_large'] = corrected_pvals_large
data['corrected_pvals_small'] = corrected_pvals_small


# add the log of corrected p-values to the dataframe
data['log_corrected_pvals_large'] = -np.log10(corrected_pvals_large)
data['log_corrected_pvals_small'] = -np.log10(corrected_pvals_small)


data.loc[data["start"] == 45769240, 'label_large'] = 'SIX5 promotor'
data.loc[data["start"] == 45770259, 'label_large'] = "CTCF peak"
data.loc[data["start"] == 45770597, 'label_large'] = "putative SIX5 enhancer"
data.loc[data["start"] == 45810188, 'label_large'] = 'RSPA6A intron 2'

data.loc[data["start"] == 45707334, 'label_small'] = "QPCTL/FBXO46 intergenic"


# save the data to a new file
data.to_csv("normal_peaks_widepeak_stats_with_pvals.txt", sep="\t", index=False)

# save a bed file with significant peaks from the large expansion
data_large = data[data['corrected_pvals_large'] < 0.05]
data_large.to_csv("largeRE_sig_peaks.bed", sep="\t", columns=['chr', 'start', 'end'], header=False, index=False)

# save a bed file with significant peaks from the small expansion
data_small = data[data['corrected_pvals_small'] < 0.05]
data_small.to_csv("smallRE_sig_peaks.bed", sep="\t", columns=['chr', 'start', 'end'], header=False, index=False)




# add a new column for labels for individual datapoints on the volcano plot
#data.loc[data["start"] == 45749503, 'label_large'] = 'Peak 1a'
#data.loc[data["start"] == 45767633, 'label_large'] = 'Peak 2a'
#data.loc[data["start"] == 45769240, 'label_large'] = 'Peak 3a, (SIX5 promotor)'
#data.loc[data["start"] == 45770259, 'label_large'] = "Peak 4a, (CTCF peak)"
#data.loc[data["start"] == 45770597, 'label_large'] = "Peak 5a"
#data.loc[data["start"] == 45781257, 'label_large'] = 'Peak 6a,2b'
#data.loc[data["start"] == 45810188, 'label_large'] = 'Peak 7a'
#data.loc[data["start"] == 45707334, 'label_small'] = 'Peak 1b'
#data.loc[data["start"] == 45781257, 'label_small'] = 'Peak 2b'

data.loc[data["start"] == 45769240, 'label_large'] = 'SIX5 promotor'
data.loc[data["start"] == 45770259, 'label_large'] = "CTCF peak"
data.loc[data["start"] == 45770597, 'label_large'] = "putative SIX5 enhancer"
data.loc[data["start"] == 45810188, 'label_large'] = 'RSPA6A intron 2'

data.loc[data["start"] == 45707334, 'label_small'] = "QPCTL/FBXO46 intergenic"




# Plot volcano plot 
fig, ax = plt.subplots()
ax.scatter(data['fraction_actuation_diff_large'], data['log_pval_large'], color='#EC2027', label='Gen II/III Expanded')
ax.scatter(data['fraction_actuation_diff_small'], data['log_pval_small'], color='#006838', label='Gen I Expanded')
plt.axhline(y=-np.log10(0.05), color='black', linestyle='--')
plt.xlabel('Fraction actuation difference')
# set x axis limits
plt.xlim(-0.65, 0.65)
plt.ylabel('-log10(p-value)')
plt.legend()

# Add labels to the points
for i in range(len(data)):
	if not pd.isnull(data['label_large'][i]):
		ax.text(data['fraction_actuation_diff_large'][i], data['log_pval_large'][i], data['label_large'][i], fontsize=8)

#for i in range(len(data)):
#	if not pd.isnull(data['label_small'][i]):
#		ax.text(data['fraction_actuation_diff_small'][i], data['log_pval_small'][i], data['label_small'][i], fontsize=8)

# Save the plot
plt.rcParams['pdf.fonttype'] = 42
plt.savefig('volcano_plot_unadj_pval.png')
plt.savefig('volcano_plot_unadj_pval.pdf', format="pdf", bbox_inches='tight')
plt.show()


# Plot volcano plot with the corrected p values
fig, ax = plt.subplots()
ax.scatter(data['fraction_actuation_diff_large'], data['log_corrected_pvals_large'], color='#EC2027', label='Gen II/III Expanded')
ax.scatter(data['fraction_actuation_diff_small'], data['log_corrected_pvals_small'], color='#006838', label='Gen I Expanded')
plt.axhline(y=-np.log10(0.05), color='black', linestyle='--')
plt.xlabel('Fraction actuation difference')
# set x axis limits
plt.xlim(-0.65, 0.65)
plt.ylabel('-log10(p-value)')
plt.legend()

# Add labels to the points
for i in range(len(data)):
	if not pd.isnull(data['label_large'][i]):
		ax.text(data['fraction_actuation_diff_large'][i], data['log_corrected_pvals_large'][i], data['label_large'][i], fontsize=8)	
	if not pd.isnull(data['label_small'][i]):
		ax.text(data['fraction_actuation_diff_small'][i], data['log_corrected_pvals_small'][i], data['label_small'][i], fontsize=8)

# Save the plot
plt.rcParams['pdf.fonttype'] = 42
plt.savefig('volcano_plot_BHpval.png')
plt.savefig('volcano_plot_BHpval.pdf', format="pdf", bbox_inches='tight')
plt.show()
