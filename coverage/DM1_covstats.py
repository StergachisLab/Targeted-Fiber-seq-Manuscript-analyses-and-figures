import pysam_coverage as cov
import numpy as np


stats_txt_path="/mmfs1/gscratch/stergachislab/bohaczuk/analysis/targeted_methods_paper/pysam_covplots/DM1_covstats.txt"
pad=0
regions=[['chr3', 129087203, 129272320],['chr4', 3006837, 3125654],['chr19', 45664901, 45825746], ['chr14', 20283833, 20402650]] # ['chr11', 5186600, 5304585],




samples=[['PS00150', 0.12], ['PS00151', 0.23], ['PS00152', 0.11], ['PS00153', 0.18], ['PS00442', 0.24], ['PS00443', 0.21], ['PS00444', 0.17]]





stats_file=open(stats_txt_path, "w")

med_covs=[]


for sample in samples:
	stats_file.write("Sample: " + sample[0] + "\n")
	bampath="/mmfs1/gscratch/stergachislab/bohaczuk/data/ft-m6a/" +  sample[0] + "/read_filtering/" + sample[0] + ".hg38.75nuc.m6afil.bam"
	print(bampath)
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

