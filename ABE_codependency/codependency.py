"""
Pipeline to take input of peak regions in BED format and Fiber-seq MSP data and produce co-accessibility and codependency data and intermediate files.

Inputs:
    Peak regions (BED format)
    Fiber-seq read coordinates (BED or BAM format)
    MSP model results (BED format)
"""

import os
import sys
import subprocess
import pandas as pd
import numpy as np
import csv
import gzip
import itertools
from datetime import datetime
import pytz
import logging
import argparse
import scipy.stats as stats
from colour import Color


# parse command line arguments
parser = argparse.ArgumentParser(description = "Pipeline to take input of peak regions in BED format and Fiber-seq MSP data and produce co-accessibility and codependency data and intermediate files.")
parser.add_argument("-m", "--model", required = True, metavar = '', help = "accessibility model MSP results (BED format)")
parser.add_argument("-p", "--peaks", required = True, metavar = '', help = "peak regions (BED format)")
parser.add_argument("-f", "--fibers", required = True, metavar = '', help = "Fiber-seq reads (BED or BAM format)")
parser.add_argument("-o", "--output", required = True, metavar = '', help = "output directory")
parser.add_argument("-n", "--name", required = False, metavar = '', help = "name to prefix on output files. Defaults to file basename.")
parser.add_argument("-g", "--genome_sizes", required = False, metavar = '', help = "Chromosome sizes of analysis genome assembly. Defaults to hg38.")
parser.add_argument("-c", "--cov", required = False, metavar = '', help = "minimum coverage for each peak pair")
args = parser.parse_args()

output_dir = args.output
model = args.model
peaks = args.peaks
fibers = args.fibers
if args.name:
    sample_name = args.name
else:
    sample_name = os.path.basename(fibers.rstrip(".bed.gz").rstrip('.bed').rstrip('.bam'))
if args.genome_sizes:
    gsizes = args.genome_sizes
else:
    gsizes = '/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes'
gsizes = os.path.abspath(gsizes)

# min fiber coverage
if args.cov:
    try:
        int(args.cov)
    except:
        sys.exit(f"Coverage cutoff (-c, --cov) must be an integer. Value provided: {args.cov}\n")
    min_cutoff = int(args.cov)
else:
    min_cutoff = 30

# verify Fiberseq file extension is BED or BAM
if fibers.endswith('.bed') or fibers.endswith('.bed.gz'):
        fibBED = True
elif fibers.endswith('.bam'):
        fibBED = False
else:
        print(f'Unrecognized file extension {os.path.splitext(fibers)[1]} for fibers file {fibers}')
        sys.exit('Aborting run......')                                                                                                                                                                                                       

# if output directory directory doesn't exist, create it
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)
# ensure directory was successfully created (failure would likely be due to lack of permissions), exit if not
if not os.path.isdir(output_dir):
    sys.exit(f'Failed to create output directory {output_dir}. Exiting...')


# format log file
logging.basicConfig(filename='{}/{}_codep_pipeline_analysis_{}.log'.format(output_dir, sample_name, datetime.now(pytz.timezone('US/Pacific')).strftime('%Y-%m-%d')),
    level=logging.INFO, format='%(levelname)s:%(asctime)s:%(message)s')
# Create log header
USER = os.environ['USER']
logging.info(f'\nCodep pipeline is being run by {USER}.\n{output_dir} will be the output directory.\n')
logging.info(f'Sample name: {sample_name}')
logging.info(f'Genome assembly sizes file: {gsizes}\n')
logging.info('The following reads input files will be used:')
logging.info(f'    Model: {model}')
logging.info(f'    Peak set: {peaks}')
logging.info(f'    Fiber-seq reads: {fibers}\n')


def call_command(shell_command):
    # call shell command and log output
    try:
        output = subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)

        # format log message
        if len(output.stdout) == 0 and len(output.stderr) == 0:
            log_command = 'Shell command {} returned code {}\n'.format(output.args, output.returncode)
        elif len(output.stdout) > 0 and len(output.stderr) == 0:
            log_command = 'Shell command {} returned code {}. \nOutput: {}\n'.format(output.args, output.returncode, output.stdout)
        elif len(output.stderr) > 0 and len(output.stdout) == 0:
            log_command = 'Shell command {} returned code {}. \nStderr: {}\n'.format(output.args, output.returncode, output.stderr)
        else:
            log_command = 'Shell command {} returned code {}. \nOutput: {}\nStderr: {}\n'.format(output.args, output.returncode, output.stdout, output.stderr)

        logging.info(log_command)
        return(0)
    except Exception as error:
        logging.exception(error)
        return(1)


# step 1: compute overlaps
fiber_overlaps_tsv = f'{output_dir}/{sample_name}_fiber_overlaps.tsv'
msp_overlaps_bed = f'{output_dir}/{sample_name}_model_intersect.bed.gz'

# per fiber total peak overlaps
print("Computing per fiber peak overlaps...")
awk_statement_bed = '\' BEGIN { OFS="\\t" } { if($NF >= 46) { print $4,$1,$8,$9,$NF } }\''
awk_statement_bam = '\' BEGIN { OFS="\\t" } { if($NF >= 46) { print $4,$13,$14,$15,$NF } }\''
if fibBED:
        command = f'bedtools intersect -wo -a {fibers} -b {peaks} | awk {awk_statement_bed} > {fiber_overlaps_tsv}'
else:
        command = f'bedtools intersect -bed -wo -a {fibers} -b {peaks} | awk {awk_statement_bam} > {fiber_overlaps_tsv}'
exit_code = call_command(command)
print("Peak overlaps complete\n")

# per fiber MSP overlaps
print("Computing per fiber MSP overlaps...")
awk_statement = '\'{ if($10 != 1.0) print }\''
command = f'bedtools intersect -wo -a {model} -b {peaks} | awk {awk_statement} | bgzip -@ 10 > {msp_overlaps_bed}'
exit_code = call_command(command)
print("MSP overlaps complete\n")


# setp 2: count per-fiber peak overlaps and classify regions as actuated or not
fdr_cutoff = 0.05

encodings = f'{output_dir}/{sample_name}_{os.path.basename(peaks).rstrip(".bed.gz")}_peak_encodings.tsv'
co_acc_outfile = f'{output_dir}/{sample_name}_co-acc_{fdr_cutoff}_cutoff.tsv'

logging.info(f'Computing pre-fiber peak overlaps with an FDR cutoff of {fdr_cutoff}...')
logging.info(f'Output files:')
logging.info(f'    Peak encodings: {encodings}')
logging.info(f'    Overlap file: {co_acc_outfile}\n')


# create encoding for peaks
new = []
counter = 1
if peaks.endswith('.gz'):
    with gzip.open(peaks, 'rt') as gp:
        for line in gp:
            if not line.startswith('#'):
                temp = line.split()[:3]
                temp.append(counter)
                new.append(tuple(temp))
                counter += 1
else:
    with open(peaks) as gp:
        for line in gp:
            if not line.startswith('#'):
                temp = line.split()[:3]
                temp.append(counter)
                new.append(tuple(temp))
                counter += 1
with open(encodings, 'w') as fw:
    tsv_writer = csv.writer(fw, delimiter = '\t')
    for l in new:
        tsv_writer.writerow(l)
enc = dict()
with open(encodings) as en:
    tsv_reader = csv.reader(en, delimiter = '\t')
    for line in tsv_reader:
        name = ':'.join(line[0:3])
        id = line[3]
        enc[name] = id


# dictionary of fibers with peak overlaps and MSPs
fiber_calls = dict()

# list of per-fiber peaks overlaping (genomic coords only)
# zmw_id, peak chromosome, peak start, peak end, fiber overlap length
with open(fiber_overlaps_tsv) as fp:
    tsv_reader = csv.reader(fp, delimiter = '\t')
    for line in tsv_reader:
        zmw = line[0]
        name = ':'.join(line[1:4])
        id = int(enc[name])
        if zmw not in fiber_calls.keys():
            fiber_calls[zmw] = [set(),set()]
        fiber_calls[zmw][0].add(id)


# TRACK SAME MSP
# list of MSP peak overlaps per-fiber
msp_by_peak = dict()
with gzip.open(msp_overlaps_bed,'rt') as fm:
    tsv_reader = csv.reader(fm, delimiter = '\t')
    for line in tsv_reader:
        zmw = line[3]
        name = ':'.join(line[11:14])
        id = int(enc[name])
        fdr = float(line[9])
        mspID = '-'.join([zmw,line[1],line[2]])
        if fdr <= fdr_cutoff:
            if zmw in fiber_calls.keys():
                fiber_calls[zmw][1].add(id)
                if id in msp_by_peak.keys():
                        msp_by_peak[id].add(mspID)
                else:
                        msp_by_peak[id] = {mspID}


# for each combination of peak overlaps, add combo to dict, and count occurances and co-accessibility
# encode MSP as [ neither, first, second, both ]
co_acc_dict = dict()
for v in fiber_calls.values():
    if len(v[0]) > 1:
        over_comb = tuple(itertools.combinations(sorted(v[0]), 2))
        acc_comb = tuple(itertools.combinations(sorted(v[1]), 2))
        for c in over_comb:
            key = f'{c[0]}:{c[1]}'
            if key not in co_acc_dict.keys():
                co_acc_dict[key] = [0,0,0,0]
            if c in acc_comb:
                co_acc_dict[key][3] += 1
            elif c[0] in v[1]:
                co_acc_dict[key][1] += 1
            elif c[1] in v[1]:
                co_acc_dict[key][2] += 1
            else:
                co_acc_dict[key][0] += 1


# write output to TSV
outl = []
for k,v in co_acc_dict.items():
    temp = v.copy()
    temp.insert(0,k)
    outl.append(temp)
with open(co_acc_outfile, 'w') as fw:
    tsv_writer = csv.writer(fw, delimiter = '\t')
    for row in outl:
        tsv_writer.writerow(row)


# step 3: compute co-accessibility & codependency scores, write results to TSV
co_acc_scores_out = f'{output_dir}/{sample_name}_scores_{fdr_cutoff}.tsv'

logging.info(f'Computing co-accessibility & codependency scores...')
logging.info(f'Output file:')
logging.info(f'    Score file: {co_acc_scores_out}\n')

counts = pd.read_csv(co_acc_outfile, sep='\t',
    header = None,
    names = ['encoding', 'neither', 'first', 'second', 'both'])

# add total coverage
counts['total'] = counts.iloc[:,1:5].sum(axis=1)

# ADD NUMBER OF SAME MSP for both peaks
same_msp = []
for i in range(len(counts)):
    e1,e2 = counts.iloc[i]['encoding'].split(':')
    try:
        num_same_msp = len(msp_by_peak[int(e1)].intersection(msp_by_peak[int(e2)]))
        same_msp.append(num_same_msp)
    except:
        same_msp.append(0)
counts['same_msp'] = same_msp

# add peak regions to counts df
enc2 = dict()
with open(encodings) as fr:
    tsv_reader = csv.reader(fr, delimiter = '\t')
    for line in tsv_reader:
        enc2[int(line[3])] = {'chr':line[0], 'start':int(line[1]), 'end':int(line[2])}

enc_split = counts['encoding'].str.split(':')
p1 = []
p2 = []
for e in enc_split:
    p1.append(int(e[0]))
    p2.append(int(e[1]))

peak_l = []
for p in p1:
    peak_l.append([enc2[p]['chr'], enc2[p]['start'], enc2[p]['end']])
temp_df = pd.DataFrame(peak_l, columns = ['p1_chr','p1_start','p1_end'])
counts = pd.concat([counts, temp_df], axis=1)
peak_l = []
for p in p2:
    peak_l.append([enc2[p]['chr'], enc2[p]['start'], enc2[p]['end']])
temp_df = pd.DataFrame(peak_l, columns = ['p2_chr','p2_start','p2_end'])
counts = pd.concat([counts, temp_df], axis=1)
# add distance between peak centers
p1_centers = counts['p1_start'] + ((counts['p1_end'] - counts['p1_start']) / 2)
p2_centers = counts['p2_start'] + ((counts['p2_end'] - counts['p2_start']) / 2)
counts['dist_cen'] = p2_centers - p1_centers
# add distance between inside peak edges
counts['dist_edges'] = counts['p2_start'] - counts['p1_end']
# add ratios of accessibility for each individual peak
counts['rat1'] = (counts['first'] + counts['both']) / counts['total']
counts['rat2'] = (counts['second'] + counts['both']) / counts['total']
# ratio of both accessible is Co-accessibility
counts['co-acc'] = counts['both'] / counts['total']
# calculate codependency (difference between observed and expected co-accessibility)
counts['expected'] = counts['rat1'] * counts['rat2']
# score for codependency
counts['codep'] = counts['co-acc'] - counts['expected']

# adjusted codependency score WITHOUT counts from MSPs overlapping both peaks
adj_total = counts['total']-counts['same_msp']
adj_both = counts['both']-counts['same_msp']
adj_rat1 = (counts['first'] + adj_both) / adj_total
adj_rat2 = (counts['second'] + adj_both) / adj_total
adj_coacc = adj_both/adj_total
adj_expected = adj_rat1 * adj_rat2
counts['adj_codep'] = adj_coacc - adj_expected

# compute Fisher Exact Test p-values
pv_l = []
for index, row in counts.iterrows():
    table = [[row['neither'],row['first']], [row['second'],row['both']]]
    odds_rat, pv = stats.fisher_exact(table, alternative='two-sided')
    pv_l.append(pv)
counts['fisher_pv'] = pv_l

# scale codependency score to between 1 & -1
counts['codep'] = counts['codep'] * 4
counts['adj_codep'] = counts['adj_codep'] * 4

# filter peak pairs by coverage
counts = counts[counts['total'] >= min_cutoff]


# write results to TSV
counts.to_csv(co_acc_scores_out,
    sep='\t',
    index = False)


# move intermediate files to a new directory
int_dir = f'{output_dir}/intermediate'
if not os.path.isdir(int_dir):
        os.mkdir(int_dir)

os.rename(fiber_overlaps_tsv, f'{int_dir}/{sample_name}_fiber_overlaps.tsv')
os.rename(msp_overlaps_bed, f'{int_dir}/{sample_name}_model_intersect.bed.gz')
os.rename(co_acc_outfile, f'{int_dir}/{sample_name}_co-acc_{fdr_cutoff}_cutoff.tsv')


# create Codependency Interact Track

# filter to +/- 0.1
cutoff = 0.1
filtered = counts[abs(counts.codep) >= cutoff]

# filter by actuation
filtered = filtered[filtered.rat1 >= 0.15]
filtered = filtered[filtered.rat2 >= 0.15]

# get encodings for names
enc_split = filtered['encoding'].str.split(':', expand=True)

# reformat columns
new = pd.DataFrame(filtered[['p1_chr','p1_start','p2_end']])
new['name'] = filtered['encoding']
new['score'] = ((filtered['codep']*500)+500).astype(int)
new['value'] = filtered['codep']
new['exp'] = '.'
new['color'] = '0'
new['sourceChrom'] = filtered['p1_chr']
new['sourceStart'] = filtered['p1_start']
new['sourceEnd'] = filtered['p1_end']
new['sourceName'] = enc_split[0]
new['sourceStrand'] = '.'
new['targetChrom'] = filtered['p2_chr']
new['targetStart'] = filtered['p2_start']
new['targetEnd'] = filtered['p2_end']
new['targetName'] = enc_split[1]
new['targetStrand'] = '.'

# assign color spectrum
gray = '#808080'
lightred = '#ffcccb'
darkred = '#ff0000'
blue = '#0000ff'
lightblue = '#99ccff'
darkblue = '#00008b'

c_gray = Color(gray)
# color gradient for positive coop scores
c1 = Color(lightred)
c2 = Color(darkred)
red_range = list(c1.range_to(c2, 450))
# negative coop scores
c1 = Color(darkblue)
c2 = Color(lightblue)
blue_range = list(c1.range_to(c2, 450))

full_grad_hex = blue_range + [c_gray for i in range(100)] + red_range
for i in range(len(full_grad_hex)):
    full_grad_hex[i] = str(full_grad_hex[i])
full_grad_hex = [s.replace('gray', gray) for s in full_grad_hex]
full_grad_hex = [s.replace('red', darkred) for s in full_grad_hex]
full_grad_hex = [s.replace('DarkBlue', darkblue) for s in full_grad_hex]

# remove incorrect length 3 hex codes generated for subset of 'colour.Color' gradient values
for i in range(len(full_grad_hex)):
    if len(str(full_grad_hex[i])) != 7:
        full_grad_hex[i] = full_grad_hex[i-1]
hex_dict = dict()
for i in range(len(full_grad_hex)):
    hex_dict[i] = full_grad_hex[i]

# update colors
new['color'] = new['score'].map(hex_dict)

# write outfile
int_name = f"{output_dir}/{sample_name}_interact_track.bed"
new.to_csv(int_name, sep = '\t', index = False, header = False)

sorted_name = int_name.replace('.bed', '_SORTED.bed')
command = f'sort-bed {int_name} > {sorted_name}'
exit_code = call_command(command)

interact_as = '/mmfs1/gscratch/stergachislab/swansoe/interact.as'
command = f"/gscratch/stergachislab/install_dir/bedToBigBed -as={interact_as} -type=bed5+13 {sorted_name} {gsizes} {int_name.replace('.bed','.bb')}"
exit_code = call_command(command)

# add head to BED to make file readable in UCSC browser
head1 = f'track type=interact name="{sample_name} Codependency Interact" description="{sample_name} Codependency Interact interact" Directional=true maxHeightPixels=300:150:50 visibility=full\n'
head2 = '#chrom  chromStart  chromEnd  name  score  value  exp  color  sourceChrom  sourceStart  sourceEnd  sourceName  sourceStrand  targetChrom  targetStart  targetEnd  targetName  targetStrand\n'

new_interact_lines = []
new_interact_lines.append(head1)
new_interact_lines.append(head2)
with open(sorted_name) as fr:
    for row in fr:
        new_interact_lines.append(row)

# delete intermediate files
os.remove(int_name)
os.remove(sorted_name)

with open(int_name, 'w') as fw:
    for row in new_interact_lines:
        fw.write(row)


print('FINSIHED')
