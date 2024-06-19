#Checks for statistically significant difference between edited and unedited alleles at each peak center

sample1="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/unedited_0.15_55nuc_fire.4/fire/unedited_0.15_55nuc_fire.4.fire.bam" # path to fire bam, output of fiberseq-fire
sample2="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/edited_0.15_55nuc_fire.4/fire/edited_0.15_55nuc_fire.4.fire.bam" 
output_prefix="edited_vs_unedited_"
chr=chr11
peak_path="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/edited_0.15_55nuc_fire.4/FDR-peaks/FDR-FIRE-peaks.bed.gz" # path to peaks bed output of fiberseq-fire


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

sample1_FIREs=$(mktemp)
sample2_FIREs=$(mktemp)


samtools view -bu $sample1 $chr:$start_peak-$end_peak | ft fire --extract > "$sample1_FIREs"
samtools view -bu $sample2 $chr:$start_peak-$end_peak | ft fire --extract > "$sample2_FIREs"

conda deactivate

conda activate fiber-views





for peak in "${peaks[@]}"; do

peak_start=$(echo "$peak" | awk '{print $1}')
peak_end=$(echo "$peak" | awk '{print $2}')

sample1_FIRE_IDs=$(mktemp)
sample1_other_IDs=$(mktemp)

sample2_FIRE_IDs=$(mktemp)
sample2_other_IDs=$(mktemp)



awk -v peak_start="$peak_start" -v peak_end="$peak_end" '$2 < peak_end && $3 > peak_start' $sample1_FIREs | awk '$5<11 {print $4}' | sort |uniq > "$sample1_FIRE_IDs"
awk -v peak_start="$peak_start" -v peak_end="$peak_end" '$2 < peak_end && $3 > peak_start' $sample1_FIREs | awk '$5>11 {print $4}'| sort | uniq > "$sample1_other_IDs"


dup_count=$(comm -12 $sample1_FIRE_IDs $sample1_other_IDs| wc -l)


sample1_FIRE_count=$(wc -l < $sample1_FIRE_IDs)
sample1_other_count_unadj=$(wc -l < "$sample1_other_IDs")
sample1_other_count=$((sample1_other_count_unadj - dup_count))
#echo $peak_start $peak_end $sample1_FIRE_count $sample1_other_count_unadj $dup_count $sample1_other_count 

awk -v peak_start="$peak_start" -v peak_end="$peak_end" '$2 < peak_end && $3 > peak_start' $sample2_FIREs | awk '$5<11 {print $4}' | sort |uniq > "$sample2_FIRE_IDs"
awk -v peak_start="$peak_start" -v peak_end="$peak_end" '$2 < peak_end && $3 > peak_start' $sample2_FIREs | awk '$5>11 {print $4}'| sort | uniq > "$sample2_other_IDs"

sample2_dup_count=$(comm -12 $sample2_FIRE_IDs $sample2_other_IDs| wc -l)

sample2_FIRE_count=$(wc -l < $sample2_FIRE_IDs)
sample2_other_count_unadj=$(wc -l < "$sample2_other_IDs")
sample2_other_count=$((sample2_other_count_unadj - sample2_dup_count))
#echo $peak_start $peak_end $sample2_FIRE_count $sample2_other_count_unadj $sample2_dup_count $sample2_other_count 




printf "$chr\t$peak\t$sample1_FIRE_count\t$sample1_other_count\t$sample2_FIRE_count\t$sample2_other_count\n" >> peak_stats.txt



done


printf 'chr\tstart\tend\tsample1_FIRE\tsample1_other\tsample2_FIRE\tsample2_other\n' > ${output_prefix}counts.txt

sort peak_stats.txt | uniq >> ${output_prefix}counts.txt




rm $sample1_FIREs
rm $sample2_FIREs

rm peak_stats.txt
