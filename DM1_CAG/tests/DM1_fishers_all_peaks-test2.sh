#Checks for statistically significant difference between normal and RE alleles at each peak center

DM1_large_RE="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/DM1all.RE.4/fire/DM1all.RE.4.fire.bam"
DM1_small_RE="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/PS00442.RE.4/fire/PS00442.RE.4.fire.bam"
DM1_normal="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/DM1all.norm.4/fire/DM1all.norm.4.fire.bam"
output_prefix="normal_peaks_wide"
chr=chr19
peak_path="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/DM1all.norm.4/FDR-peaks/FDR-FIRE-peaks.bed.gz"


readarray -t peaks < <(zcat "${peak_path}" | tail -n +2 | awk 'BEGIN{OFS="\t"} {printf "%.0f\t%.0f\n", $2, $3}')



source /gscratch/stergachislab/bohaczuk/tools/Miniconda3/etc/profile.d/conda.sh
conda activate ft-pipeline


# Compare peaks


first_line="${peaks[0]}"
start_peak=$(echo "$first_line" | awk '{print $1}')

# Extract the second coordinate of the last line
last_line="${peaks[-1]}"
end_peak=$(echo "$last_line" | awk '{print $2}')


echo $start_peak
echo $end_peak

norm_FIREs=$(mktemp)
large_FIREs=$(mktemp)
small_FIREs=$(mktemp)

samtools view -bu $DM1_normal $chr:$start_peak-$end_peak | ft fire --extract > "$norm_FIREs"
samtools view -bu $DM1_large_RE $chr:$start_peak-$end_peak | ft fire --extract > "$large_FIREs"
samtools view -bu $DM1_small_RE $chr:$start_peak-$end_peak | ft fire --extract > "$small_FIREs"

conda deactivate

conda activate fiber-views

#printf 'DM1 locus promoter actuation\n' > peak_stats.txt

printf '$chr\tpos\tnorm_FIRE\tnorm_other\tlarge_FIRE\tlarge_other\tpval_large\tsmall_FIRE\tsmall_other\tpval_small\n' >> peak_stats.txt


for peak in "${peaks[@]}"; do

peak_start=$(echo "$peak" | awk '{print $1}')
peak_end=$(echo "$peak" | awk '{print $2}')


norm_FIRE_count=$(awk -v peak_start="$peak_start" -v peak_end="$peak_end" '$2 < peak_end && $3 > peak_start' $norm_FIREs | awk '$5<11'| wc -l)
norm_other_count=$(awk -v peak_start="$peak_start" -v peak_end="$peak_end" '$2 < peak_end && $3 > peak_start' $norm_FIREs | awk '$5>11'| wc -l)

large_RE_FIRE_count=$(awk -v peak_start="$peak_start" -v peak_end="$peak_end" '$2 < peak_end && $3 > peak_start' $large_FIREs | awk '$5<11'| wc -l)
large_RE_other_count=$(awk -v peak_start="$peak_start" -v peak_end="$peak_end" '$2 < peak_end && $3 > peak_start' $large_FIREs | awk '$5>11'| wc -l)

small_RE_FIRE_count=$(awk -v peak_start="$peak_start" -v peak_end="$peak_end" '$2 < peak_end && $3 > peak_start' $small_FIREs | awk '$5<11'| wc -l)
small_RE_other_count=$(awk -v peak_start="$peak_start" -v peak_end="$peak_end" '$2 < peak_end && $3 > peak_start' $small_FIREs | awk '$5>11'| wc -l)


pval_norm_large=$(python -c "from scipy.stats import fisher_exact; fisher=fisher_exact([[$large_RE_FIRE_count, $large_RE_other_count], [$norm_FIRE_count, $norm_other_count]]); print(fisher[1])")


pval_norm_small=$(python -c "from scipy.stats import fisher_exact; fisher=fisher_exact([[$small_RE_FIRE_count, $small_RE_other_count], [$norm_FIRE_count, $norm_other_count]]); print(fisher[1])")


printf "$chr\t$peak\t$norm_FIRE_count\t$norm_other_count\t$large_RE_FIRE_count\t$large_RE_other_count\t$pval_norm_large\t$small_RE_FIRE_count\t$small_RE_other_count\t$pval_norm_small\n" >> peak_stats.txt


sort peak_stats.txt | uniq > ${output_prefix}peak_stats.txt


done


rm $norm_FIREs
rm $large_FIREs
rm $small_FIREs
#rm peak_stats.txt
