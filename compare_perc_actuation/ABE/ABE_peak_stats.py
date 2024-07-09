import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy



# Load data. Firs line is the header
data = pd.read_csv("edited_vs_unedited_counts.txt", sep="\t", header=0)



# Calculate fraction actuation for each line and save as a new column
data['fraction_actuation_unedited']= data['sample1_FIRE']/(data['sample1_FIRE']+data['sample1_other'])
data['fraction_actuation_edited']= data['sample2_FIRE']/(data['sample2_FIRE']+data['sample2_other'])


# Calculate the difference between the fraction actuation, save as a new column
data['fraction_actuation_diff'] = data['fraction_actuation_edited'] - data['fraction_actuation_unedited']

# Calculate the % induction of editing
data['percent_induction'] = (data['fraction_actuation_edited'] - data['fraction_actuation_unedited'])/data['fraction_actuation_unedited']*100

# Calculate the p-value with Fisher's exact test  and save as a new column
data['pval'] = data.apply(lambda row: scipy.stats.fisher_exact([[row['sample1_FIRE'], row['sample1_other']], [row['sample2_FIRE'], row['sample2_other']]])[1], axis=1)
print(data['pval'])

# Calculate -log10(p-value) and save as a new column
data['log_pval'] = -np.log10(data['pval'])


# apply benjamini-hochberg correction
corrected_pvals=scipy.stats.false_discovery_control(data['pval'])
print(corrected_pvals)

# add the corrected p-values to the dataframe
data['corrected_pvals'] = corrected_pvals



# add the log of corrected p-values to the dataframe
data['log_corrected_pvals'] = -np.log10(corrected_pvals)

# add labels to the HBG promoters

data.loc[data["start"] == 5249784, 'label'] = 'HBG1 promoter'
data.loc[data["start"] == 5254485, 'label'] = 'HBG2 promoter'


# save the data to a new file
data.to_csv("ABE_peak_stats_with_pvals.txt", sep="\t", index=False)

# save a bed file with significant peaks
data_sig = data[data['corrected_pvals'] < 0.05]
data_sig.to_csv("ABE_sig_peaks.bed", sep="\t", columns=['chr', 'start', 'end'], header=False, index=False)






# Plot volcano plot 
fig, ax = plt.subplots()
ax.scatter(data['fraction_actuation_diff'], data['log_pval'])
plt.axhline(y=-np.log10(0.05), color='black', linestyle='--')
ax.text(-0.6, -np.log10(0.05)+0.1, 'Nominal significance (p=0.05)', fontsize=8)
plt.xlabel('Fraction actuation difference')
# set x axis limits
plt.xlim(-0.65, 0.65)
plt.ylabel('-log10(p-value)')

# add a line for Benjamini-Hochberg corrected p-value threshold
largest_sig_index = len([x for x in data['corrected_pvals'] if x < 0.05])
print(largest_sig_index)
print(-np.log10(0.05*(largest_sig_index+1)/len(data)))
print(-np.log10(.00909))
pseudo_threshold=-np.log10(0.05*(largest_sig_index+1)/len(data))
plt.axhline(y=pseudo_threshold, color='red', linestyle='--')
ax.text(-0.6, pseudo_threshold + 0.1, 'FDR-adjusted significance', fontsize=8, color='red')


# Add labels to the points
for i in range(len(data)):
	if not pd.isnull(data['label'][i]):
		ax.text(data['fraction_actuation_diff'][i], data['log_pval'][i], data['label'][i], fontsize=8)

# Save the plot
plt.rcParams['pdf.fonttype'] = 42
plt.savefig('ABE_volcano_plot_unadj_pval.png')
plt.savefig('ABE_volcano_plot_unadj_pval.pdf', format="pdf", bbox_inches='tight')
plt.show()


# Plot volcano plot with the corrected p values
fig, ax = plt.subplots()
ax.scatter(data['fraction_actuation_diff'], data['log_corrected_pvals'])
plt.axhline(y=-np.log10(0.05), color='black', linestyle='--')
plt.xlabel('Fraction actuation difference')
# set x axis limits
plt.xlim(-0.65, 0.65)
plt.ylabel('-log10(p-value)')

# Add labels to the points
for i in range(len(data)):
	if not pd.isnull(data['label'][i]):
		ax.text(data['fraction_actuation_diff'][i], data['log_corrected_pvals'][i], data['label'][i], fontsize=8)

# Save the plot
plt.rcParams['pdf.fonttype'] = 42
plt.savefig('ABE_volcano_plot_BHpval.png')
plt.savefig('ABE_volcano_plot_BHpval.pdf', format="pdf", bbox_inches='tight')
plt.show()




# plot the % induction of editing as a bar plot for peaks 3-9, with the start position as the x-axis
fig, ax = plt.subplots()
data['percent_induction'][4:10].plot(kind='bar', color='blue')

# output for bedgraph 
data.to_csv("ABE_percent_induction.bedgraph", sep="\t", columns=['chr', 'start', 'end', 'percent_induction'], header=False, index=False)

