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

# Calculate the p-value with Fisher's exact test  and save as a new column
data['pval'] = data.apply(lambda row: scipy.stats.fisher_exact([[row['sample1_FIRE'], row['sample1_other']], [row['sample2_FIRE'], row['sample2_other']]])[1], axis=1)

# Calculate -log10(p-value) and save as a new column
data['log_pval'] = -np.log10(data['pval'])


# apply benjamini-hochberg correction
corrected_pvals=scipy.stats.false_discovery_control(data['pval'])

# add the corrected p-values to the dataframe
data['corrected_pvals'] = corrected_pvals


# add the log of corrected p-values to the dataframe
data['log_corrected_pvals'] = -np.log10(corrected_pvals)


# save the data to a new file
data.to_csv("ABE_peak_stats_with_pvals.txt", sep="\t", index=False)

# save a bed file with significant peaks
data_sig = data[data['corrected_pvals'] < 0.05]
data_sig.to_csv("ABE_sig_peaks.bed", sep="\t", columns=['chr', 'start', 'end'], header=False, index=False)






# Plot volcano plot 
fig, ax = plt.subplots()
ax.scatter(data['fraction_actuation_diff'], data['log_pval'])
plt.axhline(y=-np.log10(0.05), color='black', linestyle='--')
plt.xlabel('Fraction actuation difference')
# set x axis limits
plt.xlim(-0.65, 0.65)
plt.ylabel('-log10(p-value)')

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

# Save the plot
plt.rcParams['pdf.fonttype'] = 42
plt.savefig('ABE_volcano_plot_BHpval.png')
plt.savefig('ABE_volcano_plot_BHpval.pdf', format="pdf", bbox_inches='tight')
plt.show()
