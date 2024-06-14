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
output_path='/mmfs1/gscratch/stergachislab/bohaczuk/scripts/publication_repositories/Targeted-Fiber-seq/DM1_CAG'

GM06076_bam="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/PS00442.RE.4/fire/PS00442.RE.4.fire.bam"

expanded_bam="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/DM1all.RE.4/fire/DM1all.RE.4.fire.bam"

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
#	print('number spanning both', len(all_df[all_df['spans']=='B']))


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

	#	print('count', CAG_count)
	#	print('boundaries', CAG_boundaries)	
	#	print('CAG_seq', repeat_seq[first_CAG:last_CAG+3])


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


df=make_RE_df(GM06076_bam, chrom, repeat_start, repeat_end, "grandfather")
print("GM06076 min", min(df["repeat_len"].values), "\nGM06076 max", max(df["repeat_len"].values))
print("GM06076 all", [x/3 for x in df["repeat_len"].values])

df2=df=make_RE_df(expanded_bam, chrom, repeat_start, repeat_end, "longRE")
print("LongRE min", min(df2["repeat_len"].values), "\nlongRE max", max(df2["repeat_len"].values))
print("longRE all", [x/3 for x in df2["repeat_len"].values])

