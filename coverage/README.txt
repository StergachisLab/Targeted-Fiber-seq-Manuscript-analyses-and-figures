# Coverage and enrichment were calculated as follows:
# 1) bams were filtered for m6A content and primary alignment with filter_bams.txt
# 2) For PS00118, the script PS00118_covstats.py was used to calculate coverage and produce visualizations.
# 3) For all other samples, the script DM1_covstat.py or HBB_covstat.py were used to calculate coverage
# 4) Enrichment was calculated as median sequencing depth/(mass fraction * 10) (See supplementary Table 1). Enrichment for all targets is saved in enrichment.txt. Enrichment was plotted with plot_enrichment.py