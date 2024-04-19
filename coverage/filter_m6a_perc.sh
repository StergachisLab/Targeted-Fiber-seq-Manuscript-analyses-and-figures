# Filter reads with m6A fraction bounds specified

file=$1
m6a_low_bound=$2
m6a_upper_bound=$3
output_name=$4

header=/mmfs1/gscratch/stergachislab/bohaczuk/scripts/StephanieBohaczuk_StergachisLab/m6a_content_filter/header.txt

source /gscratch/stergachislab/bohaczuk/tools/Miniconda3/etc/profile.d/conda.sh
conda activate variantcallers

ft extract --all - ${file} | tail -n +2 \
| cut -f1-4,14,15 |  cat ${header} - \
> m6a_frac.txt

conda deactivate
conda activate fiber-views

python /mmfs1/gscratch/stergachislab/bohaczuk/scripts/StephanieBohaczuk_StergachisLab/m6a_content_filter/filter_m6A_perc.py m6a_frac.txt $m6a_low_bound $m6a_upper_bound

#conda deactivate

#conda activate fiberseq-smk

samtools view -b -N include_list.txt $file > ${output_name}


#pbindex $file

#zmwfilter --exclude exclude_list.txt $file unsorted.temp.bam

#samtools sort unsorted.temp.bam > ${output_name}


#rm unsorted.temp.bam
