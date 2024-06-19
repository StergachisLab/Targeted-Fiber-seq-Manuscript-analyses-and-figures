DM1_large_RE="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/DM1all.RE.4/fire/DM1all.RE.4.fire.bam"  # path to fire bam, output of fiberseq-fire
DM1_small_RE="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/PS00442.RE.4/fire/PS00442.RE.4.fire.bam" # path to fire bam, output of fiberseq-fire
DM1_normal="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/DM1all.norm.4/fire/DM1all.norm.4.fire.bam" # path to fire bam, output of fiberseq-fire
output_prefix=TSSs_



chr=chr19


readarray -t peaks < <(awk 'BEGIN{OFS="\t"} $3=="+" {print $4} $3=="-" {print $5}' "chr19_45664901_45825746") # chr19_45664901_45825746 contains gencode TSS peaks between these positions.
echo $peaks


source /gscratch/stergachislab/bohaczuk/tools/Miniconda3/etc/profile.d/conda.sh
conda activate ft-pipeline


# Compare peaks


start_peak=${peaks[0]}
end_peak=${peaks[-1]}

norm_FIREs=$(mktemp)
large_FIREs=$(mktemp)
small_FIREs=$(mktemp)

samtools view -bu $DM1_normal $chr:$start_peak-$end_peak | ft fire --extract > "$norm_FIREs"
samtools view -bu $DM1_large_RE $chr:$start_peak-$end_peak | ft fire --extract > "$large_FIREs"
samtools view -bu $DM1_small_RE $chr:$start_peak-$end_peak | ft fire --extract > "$small_FIREs"

conda deactivate

conda activate fiber-views


for peak in "${peaks[@]}"; do


norm_FIRE_count=$(awk -v peak="$peak" '$2 < peak && $3 > peak' $norm_FIREs | awk '$5<11'| wc -l)
norm_other_count=$(awk -v peak="$peak" '$2 < peak && $3 > peak' $norm_FIREs | awk '$5>11'| wc -l)

large_RE_FIRE_count=$(awk -v peak="$peak" '$2 < peak && $3 > peak' $large_FIREs | awk '$5<11'| wc -l)
large_RE_other_count=$(awk -v peak="$peak" '$2 < peak && $3 > peak' $large_FIREs | awk '$5>11'| wc -l)

small_RE_FIRE_count=$(awk -v peak="$peak" '$2 < peak && $3 > peak' $small_FIREs | awk '$5<11'| wc -l)
small_RE_other_count=$(awk -v peak="$peak" '$2 < peak && $3 > peak' $small_FIREs | awk '$5>11'| wc -l)



printf "$chr\t$peak\t$norm_FIRE_count\t$norm_other_count\t$large_RE_FIRE_count\t$large_RE_other_count\t$small_RE_FIRE_count\t$small_RE_other_count\n" >> peak_stats.txt

done


printf 'chr\tpos\tnorm_FIRE\tnorm_other\tlarge_FIRE\tlarge_other\tsmall_FIRE\tsmall_other\n' > ${output_prefix}counts.txt


sort peak_stats.txt | uniq >> ${output_prefix}counts.txt
 




rm $norm_FIREs
rm $large_FIREs
rm $small_FIREs
rm peak_stats.txt