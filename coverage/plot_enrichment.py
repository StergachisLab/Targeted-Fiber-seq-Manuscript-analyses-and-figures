import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

# Load the line seperated data
enrichment=open("enrichment.txt","r")
enrichment_data=enrichment.readlines()

enrichment.close()

# Remove the newline characters
enrichment_data=[float(x.strip()) for x in enrichment_data]
print(enrichment_data)

# Make a swarmplot of the data
plt.figure(figsize=(2, 4))
sns.violinplot(data=enrichment_data, showmedians=True)
sns.swarmplot(data=enrichment_data, color='black')

plt.ylabel("Enrichment over WGS", fontsize=16)
# Remove the x-axis labels
plt.xticks([])

plt.savefig("enrichment_plot2.png", bbox_inches='tight')

matplotlib.rcParams['pdf.fonttype'] = 42

plt.savefig("enrichment_plot2.pdf", format='pdf', bbox_inches='tight')