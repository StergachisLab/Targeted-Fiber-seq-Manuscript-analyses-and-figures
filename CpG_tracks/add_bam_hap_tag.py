import pysam

# bam file to add haplotype tag
bam_inp="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/DM1all.RE.4/fire/DM1all.RE.4.fire.bam"
# haplotype to set
set_hap=2
#  output bam name
bam_out="DM1_RE_haptagged.bam"


bam = pysam.AlignmentFile(bam_inp, "rb", check_sq=False)
header = bam.header.to_dict()
output_bam = pysam.AlignmentFile(bam_out, "wb", header=header)



with output_bam as out_f:
    with bam as in_f:
        for read in in_f:
            read.set_tag("HP", set_hap)
            out_f.write(read)
