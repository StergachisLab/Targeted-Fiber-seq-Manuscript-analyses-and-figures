#! /mmfs1/gscratch/stergachislab/bohaczuk/tools/Miniconda3/envs/fiber-views/bin/python3.10
# coding: utf-8

import pyft
import tqdm
import matplotlib.pyplot as plt
import numpy as np
import csv
from scipy.stats import fisher_exact


def read_bed_file_for_range(bed_file, target_chrom, target_start, target_end, offset_start=17, offset_end=30):
    """
    Read regions from a BED file, filtering for regions within a specific chromosomal range,
    and only considering a specific subregion within each BED entry.
    
    :param bed_file: Path to the BED file.
    :param target_chrom: Target chromosome.
    :param target_start: Start position of the target range.
    :param target_end: End position of the target range.
    :return: Dictionary with filtered subregions (key: region name, value: (chromosome, subregion start, subregion end, strand)).
    """
    regions = {}
    with open(bed_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) < 3:
                continue

            chrom, start, end = parts[0], int(parts[1]), int(parts[2])

            strand = parts[5] if len(parts) > 5 else '+'

            # Ensure that start is a number
            if start is None:
                raise ValueError(f"Start position is None for line: {line}")
	    # Adjust the start and end positions to only include the +17 to +30 positions

            if strand == "+":
                adjusted_start = start + offset_start
                adjusted_end = start + offset_end
            else:
                adjusted_start = end - offset_end
                adjusted_end = end - offset_start

            if chrom == target_chrom and adjusted_start >= target_start and adjusted_end <= target_end:
                region_name = f"{chrom}:{adjusted_start}-{adjusted_end}"
                regions[region_name] = (chrom, adjusted_start, adjusted_end, strand)

    return regions

def is_nucleosome(fiber, start, end):
    """
    Determine if a fiber represents a nucleosome within the specified region, allowing for a window on either side.
    
    :param fiber: Fiberdata object.
    :param start: Start position of the region.
    :param end: End position of the region.
    :return: True if the fiber is identified as a nucleosome in the region, False otherwise.
    """


    if fiber.nuc is None or fiber.nuc.reference_starts is None or fiber.nuc.reference_lengths is None:
        return False
   
    for nuc_start, nuc_length in zip(fiber.nuc.reference_starts, fiber.nuc.reference_lengths):
        if nuc_start is not None and nuc_length is not None:
            nuc_end = nuc_start + nuc_length
        # Check if nucleosome overlaps with the region of interest
        if nuc_start is not None and nuc_end >= start and nuc_start <= end:
            return True

    return False

def is_msp(fiber, start, end):
    """
    Determine if a fiber represents a MSP within the specified region, allowing for a window on either side.
    
    :param fiber: Fiberdata object.
    :param start: Start position of the region.
    :param end: End position of the region.
    :return: True if the fiber is identified as a MSP in the region, False otherwise.
    """
    if fiber.msp is None or fiber.msp.reference_starts is None or fiber.msp.reference_lengths is None:
        return False

    for msp_start, msp_length in zip(fiber.msp.reference_starts, fiber.msp.reference_lengths):
        if msp_start is not None and msp_length is not None:
            msp_end = msp_start + msp_length
        # Check if MSP overlaps with the region of interest, require the region to be fully within MSP
        if msp_start is not None and msp_end >= end and msp_start <= start:
            return True

    return False


def has_m6a_methylation(fiber, start, end):
    """
    Check if a fiber has m6a methylation marks within the specified region.
    """

    if start is None or end is None:
        raise ValueError(f"Start or end positions are None for region {start}-{end}")

    if not fiber.m6a or not fiber.m6a.reference_starts:
        return False


    # Check if any methylation marks are within the adjusted region
    for m6a_start in fiber.m6a.reference_starts:
        if m6a_start is not None and start <= m6a_start <= end:
            return True

    return False

def analyze_ctcf_binding_nucleosomes(bam_file, regions):
    """
    Analyze the binding status of CTCF in specified regions considering the DNA strand, nucleosomes, and m6a methylation.
    
    :param bam_file: Path to the BAM file.
    :param regions: Dictionary with regions of interest.
    :return: Dictionary with binding status and nucleosome count for each region.
    """
    fiberbam = pyft.Fiberbam(bam_file)
    binding_status = {region: {'bound': 0, 'unbound': 0, 'nucleosomes': 0, 'msps':0, 'uncat':0, 'fiber_count': 0, '+_strand': 0, '-_strand': 0} for region in regions}
    for region_name, (chrom, start, end, strand) in regions.items():
            fibers = fiberbam.fetch(chrom, start, end)
	   	    
            for fiber in tqdm.tqdm(fibers):
                binding_status[region_name]['fiber_count'] += 1
                
				# Increment the strand-specific counters
                if fiber.strand == '+':
                    binding_status[region_name]['+_strand'] += 1
                elif fiber.strand == '-':
                    binding_status[region_name]['-_strand'] += 1
                    
               
                if is_nucleosome(fiber, start, end):
                    binding_status[region_name]['nucleosomes'] += 1
                    
                elif is_msp(fiber, start, end):
                    binding_status[region_name]['msps'] += 1
                    if not has_m6a_methylation(fiber, start, end):
                         binding_status[region_name]['bound'] += 1
                    else:
                         binding_status[region_name]['unbound'] += 1
                else:
                	binding_status[region_name]['uncat'] += 1
    return binding_status


# Specify the target chromosomal range for filtering BED file regions
#region of binding sites to look at

target_chrom = 'chr19'
target_start = 45663977
target_end = 45826330

#input data
bed_file = 'CTG_adjacent_CTCF.bed' #bed file of binding sites
norm_bam_file = '/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/DM1all.norm.4/fire/DM1all.norm.4.fire.bam' # merged bam of all normal haplotypes, output of FIRE pipeline
re_bam_file='/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/DM1all.RE.4/fire/DM1all.RE.4.fire.bam'  # merged bam of Gen II/III expanded haplotype, output of FIRE pipeline
shortre_bam_file='/mmfs1/gscratch/stergachislab/bohaczuk/scripts/fiberseq-fire-0.0.4/FIRE/results/PS00442.RE.4/fire/PS00442.RE.4.fire.bam'  # Gen I expanded haplotype, output of FIRE pipeline


def calculate_bound(bed_file, bam_file, target_chrom, target_start, target_end):
	#get regions from given BED file
	regions = read_bed_file_for_range(
		bed_file, 
		target_chrom, 
		target_start, 
		target_end, 
		offset_start=17, # offset is added to look only at the minimal footprint for CTCF binding
		offset_end=30
	)
	# Analyze binding and nucleosomes for the filtered regions
	binding_analysis = analyze_ctcf_binding_nucleosomes(
		bam_file, 
		regions)

	# Prepare data for the table and CSV output
	headers = ["Region", "Bound Count", "Unbound Count", "Nucleosome Count", "Uncategorized Count", "Fiber Count", "+ Strand Count", "- Strand Count", "Fraction Bound", "pval"]
	table_data = [headers]
	csv_data = [headers[:-1]]


	perc_bound={}
	perc_unbound={}
	perc_nucleosomes={}
	for region, counts in binding_analysis.items():
		row = [
			region,
			counts['bound'],
			counts['unbound'],
			counts['nucleosomes'],
			counts['uncat'],
			counts['fiber_count'],
			counts['+_strand'],
			counts['-_strand']
		]
		total_categorized_count=counts['bound']+counts['unbound']+counts['nucleosomes']
		perc_bound[region]=counts['bound']/total_categorized_count
		perc_unbound[region]=counts['unbound']/total_categorized_count
		perc_nucleosomes[region]=counts['nucleosomes']/total_categorized_count
		row.append(perc_bound[region])
	

		table_data.append(row)  # Append the row to the table data
		csv_data.append(row)    # Append the row to the CSV data

	return perc_bound, perc_unbound, perc_nucleosomes, binding_analysis, csv_data

     
norm_perc_bound, norm_perc_unbound, norm_perc_nucleosomes, norm_binding, norm_csv=calculate_bound(bed_file, norm_bam_file, target_chrom, target_start, target_end)
re_perc_bound, re_perc_unbound, re_perc_nucleosomes, re_binding, re_csv=calculate_bound(bed_file, re_bam_file, target_chrom, target_start, target_end)
shortre_perc_bound, shortre_perc_unbound, shortre_perc_nucleosomes, shortre_binding, shortre_csv=calculate_bound(bed_file, shortre_bam_file, target_chrom, target_start, target_end)



# Fisher exact test for each region and print to file
with open("norm_vs_longRE_pvals.txt", "w") as f:
	# Fisher for each region vs long expansion
	for region in norm_perc_bound.keys():
		# Get the counts for the region
		norm_bound=norm_binding[region]['bound']
		norm_unbound=norm_binding[region]['unbound'] + norm_binding[region]['nucleosomes']
		re_bound=re_binding[region]['bound']
		re_unbound=re_binding[region]['unbound'] + re_binding[region]['nucleosomes']
		# Perform the Fisher exact test
		oddsratio, pvalue = fisher_exact([[norm_bound, norm_unbound], [re_bound, re_unbound]])
		f.write(f"{region} p-value: {pvalue}\n")

# Fisher for each region vs short expansion
with open("norm_vs_shortRE_pvals.txt", "w") as f:
	# Fisher for each region vs short expansion
	for region in norm_perc_bound.keys():
		# Get the counts for the region
		norm_bound=norm_binding[region]['bound']
		norm_unbound=norm_binding[region]['unbound'] + norm_binding[region]['nucleosomes']
		shortre_bound=shortre_binding[region]['bound']
		shortre_unbound=shortre_binding[region]['unbound'] + shortre_binding[region]['nucleosomes']
		# Perform the Fisher exact test
		oddsratio, pvalue = fisher_exact([[norm_bound, norm_unbound], [shortre_bound, shortre_unbound]])
		f.write(f"{region} p-value: {pvalue}\n")

# Save the table data to a CSV file
csv_file_path = 'CTCF_norm_all.csv'  
with open(csv_file_path, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(norm_csv)
    
csv_file_path = 'CTCF_RE_all.csv' 
with open(csv_file_path, 'w', newline='') as file:
	writer = csv.writer(file)
	writer.writerows(re_csv)
     
csv_file_path = 'CTCF_shortRE_all.csv'
with open(csv_file_path, 'w', newline='') as file:
	writer = csv.writer(file)
	writer.writerows(shortre_csv)



# plot regions as a bar plot sie by side
plt.rcParams['pdf.fonttype'] = 42 
ROIs=norm_perc_bound.keys()
for ROI in ROIs:
	fig=plt.figure(figsize=(4,8))
	vals=[norm_perc_bound[ROI], re_perc_bound[ROI]]
	plt.bar(x=['Normal', 'Pathogenic'], height=vals, color=['#0C0C78', '#FF45B3'])
	plt.ylim([0,0.7])
	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)
	plt.ylabel("Fraction bound", fontsize=16)
	plt.title(ROI)

	plt.savefig(ROI.replace(":","_") +'_CTCF_binding_v2.pdf', format='pdf', bbox_inches='tight')
	plt.show()

	# plot as a stacked bar plot with normal and pathogenic side by side
	fig=plt.figure(figsize=(4,8))
	vals_bound=[norm_perc_bound[ROI], re_perc_bound[ROI]]
	vals_unbound=[norm_perc_unbound[ROI], re_perc_unbound[ROI]]
	vals_nucleosomes=[norm_perc_nucleosomes[ROI], re_perc_nucleosomes[ROI]]
	plt.bar(x=['Normal', 'Pathogenic'], height=vals_bound, color='#39B54A', label='Bound')
	plt.bar(x=['Normal', 'Pathogenic'], height=vals_unbound, bottom=vals_bound, color='#2B3990', label='Unbound')
	plt.bar(x=['Normal', 'Pathogenic'], height=vals_nucleosomes, bottom=[sum(x) for x in zip(vals_bound, vals_unbound)], color='#AAAAA9', label='Nucleosomes')
    # add a legend 
	plt.legend(fontsize=16)

	plt.savefig(ROI.replace(":","_") +'_CTCF_binding_stacked_v2.pdf', format='pdf', bbox_inches='tight')


# Fisher exact test for each region and print to file
with open("norm_vs_longRE_pvals.txt", "w") as f:
	# Fisher for each region vs long expansion
	for region in norm_perc_bound.keys():
		# Get the counts for the region
		norm_bound=norm_binding[region]['bound']
		norm_unbound=norm_binding[region]['unbound'] + norm_binding[region]['nucleosomes']
		re_bound=re_binding[region]['bound']
		re_unbound=re_binding[region]['unbound'] + re_binding[region]['nucleosomes']
		# Perform the Fisher exact test
		oddsratio, pvalue = fisher_exact([[norm_bound, norm_unbound], [re_bound, re_unbound]])
		f.write(f"{region} p-value: {pvalue}\n")

# Fisher for each region vs short expansion
with open("norm_vs_shortRE_pvals.txt", "w") as f:
	# Fisher for each region vs short expansion
	for region in norm_perc_bound.keys():
		# Get the counts for the region
		norm_bound=norm_binding[region]['bound']
		norm_unbound=norm_binding[region]['unbound'] + norm_binding[region]['nucleosomes']
		shortre_bound=shortre_binding[region]['bound']
		shortre_unbound=shortre_binding[region]['unbound'] + shortre_binding[region]['nucleosomes']
		# Perform the Fisher exact test
		oddsratio, pvalue = fisher_exact([[norm_bound, norm_unbound], [shortre_bound, shortre_unbound]])
		f.write(f"{region} p-value: {pvalue}\n")
