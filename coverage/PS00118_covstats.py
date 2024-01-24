import pysam_coverage as cov
import numpy as np


bampath="/mmfs1/gscratch/stergachislab/bohaczuk/data/ft-m6a/PS00118/read_filtering/PS00118.hg38.75nuc.m6afil.bam"
nontargeted_bam="/mmfs1/gscratch/stergachislab/bohaczuk/data/ft-m6a/PS00077/m54329U_210813_020940/read_filtering/m54329U_210813_020940_75nuc.m6afil.bam"
colors=['#2178b5', '#cc4b9b', '#3a923a', '#808285']
stats_txt_path="/mmfs1/gscratch/stergachislab/bohaczuk/analysis/targeted_methods_paper/pysam_covplots/PS00118_covstats.txt"
pad=25000
CATCH_mass_fractions=[0.58, 0.58, 0.58, 0.58]
regions=[['chr4', 3006837, 3125654, 'HTT'],['chr11', 5186600, 5304585, 'HBB'],['chr19', 45664901, 45825746, 'DMPK'], ['chr14', 20283833, 20402650, 'Off_Target']]



stats_file=open(stats_txt_path, "w")

stats_file.write("Sample: PS00118\n")



# Make cov arrays and cov visualization for each targeted region
cov_arrays=[]

labels=[]

for index, region in enumerate(regions):
	left_pad_cov, coverage, right_pad_cov=cov.make_cov_array_noindel(bampath, region[0], region[1], region[2], pad)

	covfig_path="/mmfs1/gscratch/stergachislab/bohaczuk/analysis/targeted_methods_paper/pysam_covplots/PS00118_"+region[3]+"_cov2.png"

	cov_arrays.append(coverage)

	print(region[3], np.min(coverage))

	# print the lowest 20 values in the coverage array
	print(np.sort(coverage)[:20])

	#print the index of the lowest 20 values in the coverage array
	print(np.argsort(coverage)[:20])

	labels.append(region[3])

	stats_file.write(region[0]+":"+str(region[1])+"-"+str(region[2])+"\t"+region[3]+"\n"+ 
			"median_target_coverage:" + str(np.median(coverage))+"\n" + "min_target_coverage:" 
			+ str(np.min(coverage)) + "\n" + "max_target_coverage:" + str(np.max(coverage)) + "\n")

	
	#fig=cov.coverage_bar(coverage, left_pad_cov, right_pad_cov, region[0], region[1], region[2], pad, color=colors[index])
	#fig.savefig(covfig_path, format='png', bbox_inches="tight")



# Make violin plot of coverage across all targeted regions

violinfig_path="/mmfs1/gscratch/stergachislab/bohaczuk/analysis/targeted_methods_paper/pysam_covplots/PS00118_violinplot_2.svg"

fig=cov.coverage_violin(cov_arrays, labels, colors=colors)

# set range from 0-250
fig.axes[0].set_ylim([0,250])

fig.savefig(violinfig_path, format="svg", bbox_inches="tight")
fig.savefig(violinfig_path.replace(".svg", ".pdf"), format="pdf", bbox_inches="tight")


# Make a cov array for the non-targeted ref sample
ref_cov_arrays=[]

stats_file.write("Non-targeted sample: PS00077 m54329U_210813_020940\n")

for index, region in enumerate(regions):
	left_pad_cov, coverage, right_pad_cov=cov.make_cov_array_noindel(nontargeted_bam, region[0], region[1], region[2], pad)

	ref_cov_arrays.append(coverage)

	stats_file.write(region[0]+":"+str(region[1])+"-"+str(region[2])+"\t"+region[3]+"\n"+
					 "median_target_coverage:" + str(np.median(coverage))+"\n")
	

# Make relative enrichment plot

# Using average value of 10 to represent the reference sample
reference_cov=[np.array([[10]])]*len(cov_arrays)

outpath="/mmfs1/gscratch/stergachislab/bohaczuk/analysis/targeted_methods_paper/pysam_covplots/PS00118_relative_enrichment_2.svg"

fig=cov.relative_enrichment(cov_arrays, reference_cov, labels, CATCH_mass_fractions, 
						colors=colors)

# set y axis range to 0-30
fig.axes[0].set_ylim([0,30])

fig.savefig(outpath, format="svg", bbox_inches="tight")
fig.savefig(outpath.replace(".svg", ".pdf"), format="pdf", bbox_inches="tight")




