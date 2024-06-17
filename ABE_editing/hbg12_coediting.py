import pysam
import pandas as pd
import sys
from pyfaidx import Fasta
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

reference="/mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.fa"


bamfile="/mmfs1/gscratch/stergachislab/bohaczuk/analysis/hbg_lieber/ft_reruns/hg38/hg38_all2x_55_m6afil.bam"





def editing_ids (reference, bamfile, chr, start, end, expected_edit=None):
#Returns unedited, edited, indel, and off-target edited ZIDs 

	# Make dictionary with read name and edit status from reads in bamfile
	read_dict={'edited':[], 'indel':[], 'unedited':[], 'other':[]}
	bam=pysam.AlignmentFile(bamfile, 'rb')


	for column in bam.pileup (chr, start, end, truncate=True):
		refbase = reference[chr][start:end]

		for read in column.pileups:

			querybase = read.alignment.query_sequence[read.query_position]

			queryname= read.alignment.query_name


			if read.indel != 0:
				read_dict['indel'].append(queryname)

			elif querybase == expected_edit:
				read_dict['edited'].append(queryname)

			elif refbase == querybase:
				if read.alignment.query_sequence[read.query_position-4:read.query_position + 5] == reference[chr][start-4:start+5]:
					read_dict['unedited'].append(queryname)
					
	labeled_reads=set(read_dict['edited']) | set(read_dict['indel']) | set(read_dict['unedited'])

	for column in bam.pileup(chr, start-10, end+10):
		for read in column.pileups:
			queryname= read.alignment.query_name
			if queryname not in labeled_reads and queryname not in read_dict['other']:
				read_dict['other'].append(queryname)
				print('other')


	return read_dict


def fully_spanning_fibers (bamfile, chr, start, end):
	# Returns a list of read names that fully span the region of interest
	bam=pysam.AlignmentFile(bamfile, 'rb')
	fully_spanning_reads=[]
	for read in bam.fetch(chr, start, end):
		if read.reference_start <= start and read.reference_end >= end:
			fully_spanning_reads.append(read.query_name)
	return fully_spanning_reads





ref_fasta=Fasta(reference)


hbg1_edit_dict=editing_ids (ref_fasta, bamfile, "chr11", 5249969, 5249970, expected_edit='C')

hbg2_edit_dict=editing_ids (ref_fasta, bamfile, "chr11", 5254893, 5254894, expected_edit='C')

hbg1_edited=hbg1_edit_dict['edited']

hbg2_edited=hbg2_edit_dict['edited']

hbg1_unedited=hbg1_edit_dict['unedited']

hbg2_unedited=hbg2_edit_dict['unedited']

fully_spanning=fully_spanning_fibers(bamfile, "chr11", 5249969, 5254894)

hbg12_coedited=list(set(hbg1_edited) & set(hbg2_edited) & set(fully_spanning))

hbg1_only=list((set(hbg1_edited) - set(hbg2_edited)) & set(fully_spanning))

hbg2_only=list((set(hbg2_edited) - set(hbg1_edited)) & set(fully_spanning))

neither_edited=list(set(fully_spanning) & set(hbg1_unedited) & set(hbg2_unedited))



# Plot a stacked bar plot of the editing status of each fiber


data = {
	'Neither': [len(neither_edited)],
	'HBG2 only': [len(hbg2_only)],
	'HBG1 only': [len(hbg1_only)],
	'Both': [len(hbg12_coedited)]
}

df = pd.DataFrame(data)

# Plotting
ax = df.plot(kind='bar', stacked=True, color=['#43B546', '#939598', '#6D6E71', '#FCB216'], legend=False)

# Formatting
plt.axis('off') # This removes the x-tick, since there's only one category

plt.tight_layout()


plt.savefig('hbg12_coediting.svg', format="svg", bbox_inches="tight")
plt.savefig('hbg12_coediting.pdf', format="pdf", bbox_inches="tight")
plt.show()


# output a file for each category of fiber
maindir="."

with open(maindir + "/hbg12_coedited.txt", 'w') as f:
	for zid in hbg12_coedited:
		f.write(zid + '\n')

with open(maindir + "/hbg1_only.txt", 'w') as f:
	for zid in hbg1_only:
		f.write(zid + '\n')	

with open(maindir + "/hbg2_only.txt", 'w') as f:
	for zid in hbg2_only:
		f.write(zid + '\n')

with open(maindir + "/neither_edited.txt", 'w') as f:
	for zid in neither_edited:
		f.write(zid + '\n')



# Save to file stats of HBG1 percent edited, HBG2 percent edited, HBG1 percent edited among fibers spanning both, HBG2 percent edited among fibers spanning both,
# percent co-edited among fibers spanning both, percent neither edited among fibers spanning both, percent either edited among fibers spanning both, and coedit p value

hbg1_percent_edited=len(hbg1_edited)/(len(hbg1_edited) + len(hbg1_unedited))

hbg2_percent_edited=len(hbg2_edited)/(len(hbg2_edited) + len(hbg2_unedited))

hbg1_percent_edited_spanning_both=len(set(hbg1_edited) & set(fully_spanning))/len(fully_spanning)

hbg2_percent_edited_spanning_both=len(set(hbg2_edited) & set(fully_spanning))/len(fully_spanning)

coedited_percent=len(hbg12_coedited)/len(fully_spanning)

neither_percent=len(neither_edited)/len(fully_spanning)

either_percent= 1.0-neither_percent

# p value
hbg1ed_hbg2ed=len(hbg12_coedited)
hbg1ed_hbg2uned=len(hbg1_only)
hbg1uned_hbg2ed=len(hbg2_only)
hbg1uned_hbg2uned=len(neither_edited)


coedit_p=fisher_exact([[hbg1ed_hbg2ed, hbg1ed_hbg2uned], [hbg1uned_hbg2ed, hbg1uned_hbg2uned]])[1]

with open(maindir + "/hbg12_coediting_stats.txt", 'w') as f:
	f.write("HBG1 percent edited: " + str(hbg1_percent_edited) + '\n')
	f.write("HBG2 percent edited: " + str(hbg2_percent_edited) + '\n')
	f.write("HBG1 percent edited among fibers spanning both: " + str(hbg1_percent_edited_spanning_both) + '\n')
	f.write("HBG2 percent edited among fibers spanning both: " + str(hbg2_percent_edited_spanning_both) + '\n')
	f.write("Percent co-edited among fibers spanning both: " + str(coedited_percent) + '\n')
	f.write("Percent neither edited among fibers spanning both: " + str(neither_percent) + '\n')
	f.write("Percent either edited among fibers spanning both: " + str(either_percent) + '\n')
	f.write("Co-edit p value: " + str(coedit_p) + '\n')



