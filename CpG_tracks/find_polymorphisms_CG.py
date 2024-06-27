import pysam

bam_path="/mmfs1/gscratch/stergachislab/bohaczuk/analysis/targeted_methods_paper/DM1_cpg/DM1_normvsRE_happd.bam"
ref_path="/mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.fa"
output_path="exclude_polymorphisms.txt"

bam = pysam.AlignmentFile(bam_path, "rb")
ref = pysam.FastaFile(ref_path)

variant_file= open(output_path, "w")

variants = []

for column in bam.pileup("chr19", 45668882, 45823812, truncate=True, min_base_quality=0):
	
	ref_base = ref.fetch("chr19", column.pos, column.pos+1).upper()
	#ref_base_counts = 0
	base_counts= {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'indel': 0}
	#indel_counts = 0
	for read in column.pileups:	
		if read.is_del or read.is_refskip:
			base_counts['indel'] += 1
		else:
			base_counts[read.alignment.query_sequence[read.query_position]] += 1
#		elif read.alignment.query_sequence[read.query_position] == ref_base:
#			ref_base_counts += 1
#		elif read.alignment.query_sequence[read.query_position] != ref_base:
#			alt_base_counts[read.alignment.query_sequence[read.query_position]] += 1

	
	total_counts = sum(base_counts.values())

	ref_base_counts = base_counts[ref_base]	
	#total_counts = sum(alt_base_counts.values()) + indel_counts + ref_base_counts

	diff_counts = total_counts - ref_base_counts

	if total_counts == 0:
		print("No reads at position", column.pos)

	elif diff_counts/total_counts > 0.1:
		variants.append((column.pos, ref_base, base_counts))
		variant_file.write(str(column.pos-1) + "\n")
		variant_file.write(str(column.pos) + "\n")
		
		for key in base_counts.keys():
			base_counts[key] = base_counts[key]/total_counts

		if sorted(base_counts.values())[-2] < 0.1:
			print ("check this position", column.pos, ref_base,  base_counts)

		
#		ref_flanks= ref.fetch("chr11", column.pos-1, column.pos+2).upper()
#		flanks= [ref_flanks]
#		for base in base_counts:
#			if base_counts[base] > 0:
#				alt_flanks = ref_flanks[:1] + base + ref_flanks[2:]
#		if 'CG' in flanks:
#			CG_variants.append((column.pos, ref_base, base_counts))
#		else:
#			non_CG_variants.append((column.pos, ref_base, base_counts))

#print("CG_variants", CG_variants)
#print("non_CG_variants", non_CG_variants)

#print(variants)

bam.close()
ref.close()
variant_file.close()






			