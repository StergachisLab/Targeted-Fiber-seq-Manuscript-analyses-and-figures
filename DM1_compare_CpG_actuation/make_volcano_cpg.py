import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load CpG data. Firs line is the header
data_large = pd.read_csv("largeRE_haptagged.tsv", sep="\t", header=0)
data_small = pd.read_csv("smallRE_haptagged.tsv", sep="\t", header=0)


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


# Plot two graphs, one from 0-50 and one from 150 and above
plt.ylim(150, 200)
plt.savefig('cpg_volcano_plot_150_300.png')
plt.savefig('cpg_volcano_plot_150_300.pdf', format='pdf', bbox_inches='tight')
plt.ylim(0, 50)
plt.savefig('cpg_volcano_plot_0_50.png')
plt.savefig('cpg_volcano_plot_0_50.pdf', format='pdf', bbox_inches='tight')


# Show the plot
plt.ylim(0, 200)
plt.savefig('cpg_volcano_plot_all.png')
plt.savefig('cpg_volcano_plot_all.pdf', format='pdf', bbox_inches='tight')
plt.show()





# Load long RE actuation data. First line is the header
data = pd.read_csv("../compare_perc_actuation/DM1/peaks/normal_peaks_widepeak_stats.txt", sep="\t", header=0)



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

data_large.loc[data_large["start"] == 45769240, 'label_large'] = 'SIX5 promotor'
data_large.loc[data_large["start"] == 45770259, 'label_large'] = "CTCF peak"
data_large.loc[data_large["start"] == 45770597, 'label_large'] = "putative SIX5 enhancer"



data_large['fraction_actuation_diff_large'] = -data['fraction_actuation_normal'] + data['fraction_actuation_large']


# plot cpg vs actuation
fig, ax = plt.subplots()
ax.scatter(data_large['fraction_actuation_diff_large'], data_large['mean_meth_delta'], color='#EC2027', label='Gen II/III Expanded')
plt.ylabel('mCpG difference')
plt.xlabel('Fraction actuation difference')
# Add labels to the points
for i in range(len(data_large)):
	if not pd.isnull(data_large['label_large'][i]):
		ax.text(data_large['fraction_actuation_diff_large'][i], data_large['mean_meth_delta'][i], data_large['label_large'][i], fontsize=8)

plt.show()
plt.savefig('actuation_vs_cpg_large.png')
plt.savefig('actuation_vs_cpg_large.pdf', format='pdf', bbox_inches='tight')