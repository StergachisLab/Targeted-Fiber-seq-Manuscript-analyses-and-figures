import pandas as pd
import numpy as np
import sys
import pyft
import matplotlib.pyplot as plt
import seaborn as sns
import statistics

chrom="chr19"
repeat_start=45770202
repeat_end=45770267
folder='/mmfs1/gscratch/stergachislab/bohaczuk/data/ft-m6a/'
samples=['PS00150', 'PS00151', 'PS00152', 'PS00153', 'PS00442', 'PS00443', 'PS00444']
labels=['PS00150', 'PS00151', 'PS00152', 'PS00153', 'PS00442', 'PS00443', 'PS00444']
#labels=['Normal', 'Large_pathogenic', 'Small_pathogenic']
output_path='/mmfs1/gscratch/stergachislab/bohaczuk/scripts/publication_repositories/Targeted-Fiber-seq/DM1_CAG'



def make_RE_df(bampath, chrom, repeat_start, repeat_end, label="NONE"):

	fiberbam = pyft.Fiberbam(bampath)

	all_df=pd.DataFrame(columns=['rname', 'fiber_length', 'sequence', 'bampath','right_pos', 'left_pos', 'spans', 'repeat_len', 'repeat_num', 'label']) 
	for fiber in fiberbam.fetch(chrom, start=repeat_start, end=repeat_end):
		if fiber.qname not in all_df['rname'].values:
			rname=fiber.qname
			sequence=fiber.seq
			length=int(len(sequence))
			new_df=pd.DataFrame([[rname, length, sequence, bampath, label]], columns=['rname', 'fiber_length', 'sequence', 'bampath','label'])
							
			all_df=pd.concat([all_df, new_df], ignore_index=True)
		else:
			print('duplicate fiber', fiber.qname)
			continue



	# Check all reads anchored at the left of the repeat

	start_center=fiberbam.center(chrom, start=repeat_start-1, end=repeat_start, strand="+")

	for i in start_center:
		fiber=i.qname
		left_pos=i.lift_query_positions([repeat_start])[0]
		all_df.loc[:,'left_pos'][all_df['rname']==fiber]=left_pos

	# Check reads anchored at the right of the repeat

	end_center=fiberbam.center(chrom, start=repeat_end-1, end=repeat_end, strand="+")

	for i in end_center:
		key=i.qname
		right_pos=i.lift_query_positions([repeat_end])[0]
		all_df.loc[:,'right_pos'][all_df['rname']==key]=right_pos



	# Check if fiber is anchored left, right, or both
	for i in range(0,len(all_df)):
		rows_to_drop=[]
		# if fiber is anchored on both sides, neither left_pos nor right_pos will be NaN
		if all_df['left_pos'][i] == all_df['left_pos'][i]  and all_df['right_pos'][i] == all_df['right_pos'][i]:
			all_df['spans'][i]='B'
			repeat_len=all_df['right_pos'][i]-all_df['left_pos'][i]
			all_df['repeat_len'][i]=int(all_df['right_pos'][i]-all_df['left_pos'][i])
		elif np.isnan(all_df['left_pos'][i]) == False and np.isnan(all_df['right_pos'][i]) == True:
			all_df['spans'][i]='L'
			all_df['right_pos'][i]=int(all_df['fiber_length'][i])
			repeat_len=all_df['fiber_length'][i]-all_df['left_pos'][i]
			#print('name', all_df['fiber'][i], 'repeat_lenB ', repeat_len, 'right_pos ', all_df['right_pos'][i], 'left_pos ', all_df['left_pos'][i])
			all_df['repeat_len'][i]=all_df['fiber_length'][i]-all_df['left_pos'][i]
		elif np.isnan(all_df['right_pos'][i]) == False and np.isnan(all_df['left_pos'][i]) == True:
			all_df['spans'][i]='R'
			all_df['left_pos'][i]=int(0)
			all_df['repeat_len'][i]=all_df['right_pos'][i]
		else:
			rows_to_drop.append(i)
			print('dropping fiber', all_df['rname'][i])


	all_df.drop(rows_to_drop, inplace=True)
	all_df.reset_index(drop=True, inplace=True)
		

	# Get names spanning both, check for unique names

	spanning_both=all_df[all_df['spans']=='B']
	print('number spanning both', len(all_df[all_df['spans']=='B']))


	# Count CAGs

	for fiber in spanning_both['rname']:
		left_pos=spanning_both[spanning_both['rname']==fiber]['left_pos'].values[0]
		right_pos=spanning_both[spanning_both['rname']==fiber]['right_pos'].values[0]
		repeat_seq=spanning_both[spanning_both['rname']==fiber]['sequence'].values[0][left_pos:right_pos]

		# find position of first and last CAG
		first_CAG=repeat_seq.find('CAG')
		last_CAG=repeat_seq.rfind('CAG')
		CAG_stretch_len=(last_CAG+3-first_CAG)
		CAG_boundaries=round(CAG_stretch_len/3)
		all_df.loc[all_df['rname'] == fiber, 'repeat num'] = CAG_boundaries


		# Count the number of CAGs
		CAG_count=repeat_seq.count('CAG')

		print('count', CAG_count)
		print('boundaries', CAG_boundaries)	
		print('CAG_seq', repeat_seq[first_CAG:last_CAG+3])


	# 	As a sanity check, see that the number of CAGs is close to the number of CAGs between the first and last CAG
	#	and that this number is close to the number of bps between the genomic anchoring positions 
		#print('CAG_count', CAG_count, 'CAG_boundaries', CAG_boundaries)
		#print("repeat_len", (spanning_both[spanning_both['rname']==fiber]['repeat_len'].values[0]-5)/3)

	# Check for any large interruptions and poor alignment near ends

		if CAG_count-CAG_boundaries > 50 or CAG_count-CAG_boundaries < -50:
			print("Check fiber " + fiber + " for large interruptions")

		if repeat_seq[first_CAG:first_CAG+9] != 'CAGCAGCAG' or repeat_seq[last_CAG-6:last_CAG+3] != 'CAGCAGCAG':
	#		print(repeat_seq[CAGs[0]:CAGs[0]+10], repeat_seq[CAGs[-1]:CAGs[-1]+10])
			print("Check fiber " + fiber + " in sample " + label + " for poor anchoring")
			print(repeat_seq[first_CAG:first_CAG+9], repeat_seq[last_CAG-6:last_CAG+3])

	return all_df


def categorize_and_extract(df, sizes, cat_name):
	'''
	Function to categorize reads based on repeat length
	'''
	df_mask = (df['repeat num']>=sizes[0]) & (df['repeat num'] <= sizes[1])
	df.loc[df_mask, 'designation']=cat_name
	print(len(df.loc[df['designation']==cat_name]))



all_samples_df=pd.DataFrame(columns=['rname', 'fiber_length', 'sequence', 'bampath','right_pos', 'left_pos', 'spans', 'repeat_len', 'repeat num', 'label', 'RE_designation'])
for i, sample in enumerate(samples):
	bampath=folder+sample+'/read_filtering/'+sample + '.hg38.75nuc.m6afil.keepsupplement.bam'
	df=make_RE_df(bampath, chrom, repeat_start, repeat_end, labels[i])
	all_samples_df=pd.concat([all_samples_df, df], ignore_index=True)

print(all_samples_df['spans'])
all_df_spans = all_samples_df[all_samples_df['spans'] == 'B']

# Categorize reads based on repeat length
categorize_and_extract(all_df_spans, [0, 34], 'normal')
categorize_and_extract(all_df_spans, [35, 300], 'small_pathogenic')
categorize_and_extract(all_df_spans, [301, 100000000000], 'large_pathogenic')
#categorize_and_extract(all_df_spans, [35, 10000000000], 'pathogenic')



# Change labels from sample ID to individual donor ID
all_df_spans['label'] = all_df_spans['label'].replace({'PS00150': 'GM04608', 'PS00151': 'GM04608', 'PS00152': 'GM04601', 'PS00153': 'GM04602', 'PS00442': 'GM06076', 'PS00443': 'GM04608', 'PS00444': 'GM04608'})

# Print statistics for each sample and category to a file
with open(output_path + "/RE_stats.txt", 'w') as f:
	for sample in all_df_spans['label'].unique():
		for category in ['normal', 'small_pathogenic', 'large_pathogenic']:
#		for category in ['normal', 'pathogenic']:
			f.write(sample + ' ' + category + ' ' + str(len(all_df_spans[(all_df_spans['label']==sample) & (all_df_spans['designation']==category)])) + '\n')
			repeat_numbers=all_df_spans[(all_df_spans['label']==sample) & (all_df_spans['designation']==category)]['repeat num']
			if len(repeat_numbers) > 0:
				f.write('Median repeat length: ' + str(statistics.median(repeat_numbers)) + '\n')
				f.write('Minimum repeat length: ' + str(min(repeat_numbers)) + '\n')
				f.write('Maximum repeat length: ' + str(max(repeat_numbers)) + '\n')
			else:
				f.write('No reads in this category\n')
		f.write('\n\n')



# Remove categories with no data points
data_filtered = all_df_spans.groupby(['label', 'designation']).filter(lambda x: len(x) > 0)


# Plot violin plots by category of repeat length
fig, ax = plt.subplots()
#sns.violinplot(data=data_filtered, x='label', y='repeat num', hue="designation", split=False, showmeans=False, palette='Set2', showmedians=True, edgecolor='black', linewidth=1, cut=0, alpha=0.2, scale='width')
sns.swarmplot( data=data_filtered, x='label', y='repeat num', hue="designation", size=2)
plt.ylabel("CAG Repeat Number", fontsize=16)
plt.yticks(fontsize=16)

#Make y axis log scale
plt.yscale('log')
# set xtick labels
#plt.xticks([0,1,2], ['Normal\n(all samples)', 'Pathogenic\n(mother/daughters)', 'Pathogenic\n(grandfather)'])
# decrease x tick label font size
plt.xticks(fontsize=16)
# rotate x tick labels
plt.xticks(rotation=90)
#Remove x axis label
plt.xlabel('')
#plt.xticks([0,1], ['Normal', 'Pathogenic'])
plt.savefig(output_path + "/RE_violin_sample.pdf", bbox_inches='tight')
plt.savefig(output_path + "/RE_violin_sample.png", bbox_inches='tight')


# Plot violin plots by category of repeat length
fig, ax = plt.subplots()
sns.violinplot(data=all_df_spans, x='designation', y='repeat num', showmeans=False, palette='Set2', showmedians=True, edgecolor='black', linewidth=1, cut=0, alpha=0.2, scale='width')
sns.swarmplot(data=all_df_spans, x='designation', y='repeat num', color='black', size=2)
plt.ylabel("CAG Repeat Number", fontsize=16)
plt.xlabel('')
plt.yticks(fontsize=16)

#Make y axis log scale
plt.yscale('log')
plt.savefig(output_path + "/RE_violin_category.pdf", bbox_inches='tight')
plt.savefig(output_path + "/RE_violin_category.png", bbox_inches='tight')


# Still need to categorize by individual donor, separate pathogenic and non-pathogenic for each donor

