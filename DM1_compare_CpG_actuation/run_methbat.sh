# Make file with peaks
#cut -f1-3 /mmfs1/gscratch/stergachislab/bohaczuk/scripts/publication_repositories/Targeted-Fiber-seq/DM1_CAG/tests/normal_peaks_widepeak_stats.txt > peaks.txt

# Run methbat
/mmfs1/gscratch/stergachislab/bohaczuk/analysis/hbg_lieber/ft_reruns/cpg/methbat/methbat-v0.9.0-x86_64-unknown-linux-gnu/methbat profile \
	--input-prefix /mmfs1/gscratch/stergachislab/bohaczuk/analysis/targeted_methods_paper/DM1_cpg/haptagged_cpg_track \  # path to pbcpg output
	--input-regions peaks.txt \
	--output-region-profile largeRE_haptagged.tsv

/mmfs1/gscratch/stergachislab/bohaczuk/analysis/hbg_lieber/ft_reruns/cpg/methbat/methbat-v0.9.0-x86_64-unknown-linux-gnu/methbat profile \
	--input-prefix /mmfs1/gscratch/stergachislab/bohaczuk/analysis/targeted_methods_paper/DM1_cpg/PS00442/haptagged_cpg_track \
	--input-regions peaks.txt \
	--output-region-profile smallRE_haptagged.tsv
