import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load data. Firs line is the header
data = pd.read_csv("/mmfs1/gscratch/stergachislab/bohaczuk/scripts/publication_repositories/Targeted-Fiber-seq/DM1_CAG/tests/normal_peaks_widepeak_stats.txt", sep="\t", header=0)



# Calculate fraction actuation for each line and save as a new column
data['fraction_actuation_normal']= data['norm_FIRE']/(data['norm_FIRE']+data['norm_other'])
data['fraction_actuation_large'] = data['large_FIRE']/(data['large_FIRE']+data['large_other'])
data['fraction_actuation_small'] = data['small_FIRE']/(data['small_FIRE']+data['small_other'])


# Calculate the difference between the fraction actuation, save as a new column
data['fraction_actuation_diff_large'] = -data['fraction_actuation_normal'] + data['fraction_actuation_large']
data['fraction_actuation_diff_small'] = -data['fraction_actuation_normal'] + data['fraction_actuation_small']


# Calculate -log10(p-value) and save as a new column
data['log_pval_large'] = -np.log10(data['pval_large'])
data['log_pval_small'] = -np.log10(data['pval_small'])

# add a new column for labels for individual datapoints on the volcano plot
#data.loc[data["pos_start"] == 45749503, 'label_large'] = 'Peak 1a'
#data.loc[data["pos_start"] == 45767633, 'label_large'] = 'Peak 2a'
#data.loc[data["pos_start"] == 45769240, 'label_large'] = 'Peak 3a, (SIX5 promotor)'
#data.loc[data["pos_start"] == 45770259, 'label_large'] = "Peak 4a, (CTCF peak)"
#data.loc[data["pos_start"] == 45770597, 'label_large'] = "Peak 5a"
#data.loc[data["pos_start"] == 45781257, 'label_large'] = 'Peak 6a,2b'
#data.loc[data["pos_start"] == 45810188, 'label_large'] = 'Peak 7a'
#data.loc[data["pos_start"] == 45707334, 'label_small'] = 'Peak 1b'
#data.loc[data["pos_start"] == 45781257, 'label_small'] = 'Peak 2b'

data.loc[data["pos_start"] == 45769240, 'label_large'] = 'SIX5 promotor'
data.loc[data["pos_start"] == 45770259, 'label_large'] = "CTCF peak"
data.loc[data["pos_start"] == 45770597, 'label_large'] = "putative SIX5 enhancer"




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
plt.savefig('volcano_plot.png')
plt.show()