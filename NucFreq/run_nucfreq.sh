script_path="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/NucFreq"


hg38_bam="/mmfs1/gscratch/stergachislab/bohaczuk/analysis/hbg_lieber/ft_reruns/asm/for_nucfreq/alldonor1.hg38.hbb.bam" # bam with all donor 1 reads aligned to hg38
hg38_output="all.hg38.nucfreq.png"
regions="chr11:5180424-5305559"

asm_bam="/mmfs1/gscratch/stergachislab/bohaczuk/analysis/hbg_lieber/ft_reruns/asm/asm_allABE_55.bam" # bam with all donor 1 reads aligned to hifiasm assembly
asm_output="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/publication_repositories/Targeted-Fiber-seq/NucFreq/all.asm.nucfreq.png"




source /gscratch/stergachislab/bohaczuk/tools/Miniconda3/etc/profile.d/conda.sh
conda activate fiber-views # conda environment with python

python $script_path/NucPlot.py $hg38_bam $hg38_output --regions $regions

python $script_path/NucPlot.py $asm_bam $asm_output