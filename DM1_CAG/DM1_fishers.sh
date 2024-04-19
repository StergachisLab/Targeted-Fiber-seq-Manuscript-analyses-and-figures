DM1_large_RE="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/DM1all.RE.4/fire/DM1all.RE.4.fire.bam"
DM1_small_RE="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/PS00442.RE.4/fire/PS00442.RE.4.fire.bam"
DM1_normal="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/DM1all.norm.4/fire/DM1all.norm.4.fire.bam"



chr=chr19
peaks=(45769391 45770364 45770866)



source /gscratch/stergachislab/bohaczuk/tools/Miniconda3/etc/profile.d/conda.sh
conda activate ft-pipeline

echo "DM1 peak stats" > peak_stats.txt

# Compare peaks

for peak in "${peaks[@]}"; do

conda deactivate
conda activate ft-pipeline



norm_FIRE_count=$(samtools view -u $DM1_normal $chr:$peak| ft fire --extract -  | awk -v peak="$peak" '$2 < peak && $3 > peak' | awk '$5<11'| wc -l)
norm_other_count=$(samtools view -u $DM1_normal $chr:$peak| ft fire --extract -  | awk -v peak="$peak" '$2 < peak && $3 > peak' | awk '$5>11'| wc -l)


large_RE_FIRE_count=$(samtools view -u $DM1_large_RE $chr:$peak| ft fire --extract -  | awk -v peak="$peak" '$2 < peak && $3 > peak' | awk '$5<11'| wc -l)
large_RE_other_count=$(samtools view -u $DM1_large_RE $chr:$peak| ft fire --extract -  | awk -v peak="$peak" '$2 < peak && $3 > peak' | awk '$5>11'| wc -l)

small_RE_FIRE_count=$(samtools view -u $DM1_small_RE $chr:$peak| ft fire --extract -  | awk -v peak="$peak" '$2 < peak && $3 > peak' | awk '$5<11'| wc -l)
small_RE_other_count=$(samtools view -u $DM1_small_RE $chr:$peak| ft fire --extract -  | awk -v peak="$peak" '$2 < peak && $3 > peak' | awk '$5>11'| wc -l)



echo "normal FIRE count $norm_FIRE_count" >> peak_stats.txt
echo "normal_other_count $norm_other_count" >> peak_stats.txt

echo "large RE FIRE count $large_RE_FIRE_count" >>peak_stats.txt
echo "large RE other count $large_RE_other_count" >>peak_stats.txt

echo "small RE FIRE count $small_RE_FIRE_count" >>peak_stats.txt
echo "small RE other count $small_RE_other_count" >>peak_stats.txt



conda deactivate
conda activate fiber-views


pval=$(python -c "from scipy.stats import fisher_exact; fisher=fisher_exact([[$large_RE_FIRE_count, $large_RE_other_count], [$norm_FIRE_count, $norm_other_count]]); print(fisher)")
echo "large RE vs normal $pval" >> peak_stats.txt


pval=$(python -c "from scipy.stats import fisher_exact; fisher=fisher_exact([[$small_RE_FIRE_count, $small_RE_other_count], [$norm_FIRE_count, $norm_other_count]]); print(fisher)")
echo "small RE vs normal $pval" >> peak_stats.txt



done