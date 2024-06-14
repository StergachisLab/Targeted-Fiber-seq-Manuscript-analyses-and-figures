import pysam_coverage as cov
import numpy as np


stats_txt_path="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/publication_repositories/Targeted-Fiber-seq/coverage/DM1_phased_covstats.txt"
pad=0
regions=[['chr19', 45664901, 45825746]] # ['chr11', 5186600, 5304585],


samples=["/mmfs1/gscratch/stergachislab/bohaczuk/analysis/targeted_methods_paper/bam_for_FIRE/DM1/DM1all.norm.nohap.m6afil.bam",
		"/mmfs1/gscratch/stergachislab/bohaczuk/analysis/targeted_methods_paper/bam_for_FIRE/DM1/DM1all.RE.nohap.m6afil.bam",
		"/mmfs1/gscratch/stergachislab/bohaczuk/analysis/targeted_methods_paper/bam_for_FIRE/DM1/PS00442.RE.m6afil.bam"]





stats_file=open(stats_txt_path, "w")

med_covs=[]


for sample in samples:

	bampath=sample
	print(bampath)
	stats_file.write("Sample: " + sample + "\n")
	for index, region in enumerate(regions):
		left_pad_cov, coverage, right_pad_cov=cov.make_cov_array_noindel(bampath, region[0], region[1], region[2], pad)
		
		print('l', left_pad_cov, 'coverage', coverage)
		med_covs.append(str(np.median(coverage)))

		
		stats_file.write(region[0]+":"+str(region[1])+"-"+str(region[2])+"\n"+ 
			"median_target_coverage:" + str(np.median(coverage))+"\n" + "min_target_coverage:" 
			+ str(np.min(coverage)) + "\n" + "max_target_coverage:" + str(np.max(coverage)) + "\n\n")
	

for i in range(0,len(med_covs)):
	if (i+1)%4 !=0:
		stats_file.write(med_covs[i] + '\n')

