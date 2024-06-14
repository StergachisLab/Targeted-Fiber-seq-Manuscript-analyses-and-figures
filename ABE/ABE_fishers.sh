edited="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/edited_0.15_55nuc_fire.4/fire/edited_0.15_55nuc_fire.4.fire.bam"
unedited="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/unedited_0.15_55nuc_fire.4/fire/unedited_0.15_55nuc_fire.4.fire.bam"




chr=chr11
peaks=(5247698 5249346 5249947 5252404 5254157 5254536)



source /gscratch/stergachislab/bohaczuk/tools/Miniconda3/etc/profile.d/conda.sh
conda activate ft-pipeline

echo "ABE peak stats" > peak_stats.txt

# Compare peaks

for peak in "${peaks[@]}"; do

conda deactivate
conda activate ft-pipeline



unedited_FIRE_count=$(samtools view -u $unedited $chr:$peak| ft fire --extract -  | awk -v peak="$peak" '$2 < peak && $3 > peak' | awk '$5<11'| wc -l)
unedited_other_count=$(samtools view -u $unedited $chr:$peak| ft fire --extract -  | awk -v peak="$peak" '$2 < peak && $3 > peak' | awk '$5>11'| wc -l)


edited_FIRE_count=$(samtools view -u $edited $chr:$peak| ft fire --extract -  | awk -v peak="$peak" '$2 < peak && $3 > peak' | awk '$5<11'| wc -l)
edited_other_count=$(samtools view -u $edited $chr:$peak| ft fire --extract -  | awk -v peak="$peak" '$2 < peak && $3 > peak' | awk '$5>11'| wc -l)

echo peak:$peak >> peak_stats.txt

echo "unedited FIRE count $unedited_FIRE_count" >> peak_stats.txt
echo "unedited other count $unedited_other_count" >> peak_stats.txt

echo "edited FIRE count $edited_FIRE_count" >>peak_stats.txt
echo "edited other count $edited_other_count" >>peak_stats.txt



conda deactivate
conda activate fiber-views


pval=$(python -c "from scipy.stats import fisher_exact; fisher=fisher_exact([[$edited_FIRE_count, $edited_other_count], [$unedited_FIRE_count, $unedited_other_count]]); print(fisher)")
echo $pval >> peak_stats.txt


done