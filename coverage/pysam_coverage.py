import pysam
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import sys





def make_cov_array(bampath, contig, start, end, pad):
	'''Returns left pad coverage, coverage array, and right pad coverage for a given region
	For coverage bar plots, combine the three arrays into one array. For statistics of 
	coverage across the target, use just the coverage array
	
	Parameters
	bampath: path to bam file
	contig: contig name
	start: start position of the central coverage array
	end: end position of the central coverage array
	pad: number of bases to return in the rigth and left pad arrays'''

	bamfile=pysam.AlignmentFile(bampath, "rb")
	left_pad_cov=np.sum(np.array(bamfile.count_coverage(contig, start - pad, start, quality_threshold=0)), axis=0)
	coverage=np.sum(np.array(bamfile.count_coverage(contig, start, end, quality_threshold=0)), axis=0)
	right_pad_cov=np.sum(np.array(bamfile.count_coverage(contig, end, end + pad, quality_threshold=0)), axis=0)
	
	return left_pad_cov, coverage, right_pad_cov


def make_cov_array_noindel(bampath, contig, start, end, pad):
	bamfile=pysam.AlignmentFile(bampath, "rb")

	padded_start=start-pad
	padded_end=end+pad

	if padded_start < 0:
		return ValueError("Padded start position is less than 0")

	# Open the BAM file for reading.
	coverage = np.zeros(padded_end - padded_start, dtype=int)

	# Iterate over each read aligned to the specified subregion of the reference.
	for read in bamfile.fetch(contig, padded_start, padded_end):
	# Update the coverage array based on the alignment of the read to the subregion.
		coverage_slice = slice(max(read.reference_start - padded_start, 0), 
						min(read.reference_end - padded_start, padded_end - padded_start))
		coverage[coverage_slice] += 1

	left_pad_cov=coverage[:pad]
	cov=coverage[pad:-pad]
	right_pad_cov=coverage[-pad:]

	return left_pad_cov, cov, right_pad_cov



def add_scale_bar(ax, length, position=None, text=None, orientation='horizontal', **kwargs):
    """Adds a scale bar to a plot. This code from GPT4
    
    Parameters:
    - ax: The axis to which to add the scale bar.
    - length: The length of the scale bar in axis units.
    - position: The position where the scale bar starts. If None, it defaults to the lower left corner for horizontal
      and the lower right corner for vertical orientation.
    - text: The text label for the scale bar. If None, it defaults to the same as length.
    - orientation: 'horizontal' or 'vertical'
    - kwargs: Additional arguments to pass to the Line2D constructor.
    
    Returns:
    - The scale bar (as a Line2D object).
    """
    if position is None:
        if orientation == 'horizontal':
            position = [ax.get_xlim()[0], ax.get_ylim()[0]]
        else:
            position = [ax.get_xlim()[1], ax.get_ylim()[0]]

    if text is None:
        text = str(length)

    if orientation == 'horizontal':
        bar = ax.plot([position[0], position[0] + length], [position[1], position[1]], **kwargs)
        ax.text(position[0] + length / 2, position[1], text, va='bottom', ha='center')
    else:
        bar = ax.plot([position[0], position[0]], [position[1], position[1] + length], **kwargs)
        ax.text(position[0], position[1] + length / 2, text, va='center', ha='right')

    return bar



def coverage_bar(cov_array, left_pad_cov, right_pad_cov, contig, start, end, pad, outpath=None, color=None):
	'''Makes a plot of coverage across the targeted region
	Parameters
	cov_array: coverage array across the targeted region
	left_pad: coverage array across the left pad
	right_pad: coverage array across the right pad
	contig: contig name
	start: start position of the central coverage array
	end: end position of the central coverage array
	pad: number of bases in the rigth and left pad arrays'''

	padded_cov_array=np.concatenate((left_pad_cov, cov_array, right_pad_cov))

	matplotlib.rcParams['pdf.fonttype'] = 42

	if color is None:
		color="blue"
	
	fig=plt.figure(figsize=(12, 1))
	plt.ylabel("Coverage")
	plt.bar(np.arange(len(padded_cov_array)), padded_cov_array, width=1, color=color)
	plt.xticks([])
	plt.title(contig + ":" + str(start) + "-" + str(end) + " (pad=" + str(pad) + ")")


	ax = plt.gca() # Get the current axis instance
	for spine in ax.spines.values():
		spine.set_visible(False)

	x_max, y_max = ax.get_xlim()[1], ax.get_ylim()[1]
	bar_length = 10000  # adjust this as necessary
	position = [x_max - bar_length, y_max]
	add_scale_bar(ax, bar_length, position=position, text=str(int(bar_length/1000))+"kb", color='black', orientation='horizontal')

	if outpath is not None:
		plt.savefig(outpath, format="svg", bbox_inches="tight")
		plt.savefig(outpath.replace(".svg", ".pdf"), format="pdf", bbox_inches="tight")

	else:
		return fig


def coverage_violin(cov_arrays, labels, outpath=None, colors=None):
	'''Makes a violin plot of coverage across the targeted region for multiple samples
	Parameters
	cov_arrays: list of coverage arrays for each sample
	labels: list of sample names
	outpath: path to save the plot
	color: color of the violin plot (default is blue)'''

	matplotlib.rcParams['pdf.fonttype'] = 42
	fig=plt.figure(figsize=(4,4))

	if colors is None:
		colors=["blue"]*len(cov_arrays)
	#for i in range(0,len(cov_arrays)):
	#	cov_arrays[i]=cov_arrays[i][pads[i]:-pads[i]]
	ax=sns.violinplot(cov_arrays, showmeans=False, scale='width', showmedians=True, palette=colors, edgecolor="black", linewidth=1)
	ax.set_xticklabels(labels)
	plt.ylabel("Coverage")
	plt.show()

	if outpath is not None:
		plt.savefig(outpath, format="svg", bbox_inches="tight")
		plt.savefig(outpath.replace(".svg", ".pdf"), format="pdf", bbox_inches="tight")

	else:
		return fig


def relative_enrichment(cov_arrays, cov_ref_arrays, labels, CATCH_mass_fractions, outpath=None, colors=None):

	'''Makes a bar plot of relative enrichment of multiple sample compared to a reference sample.
	Each sample will be a bar with the value determined by median coverage of the sample divided by
	the median coverage of the reference sample
	
	Parameters
	'''
	matplotlib.rcParams['pdf.fonttype'] = 42
	fig=plt.figure(figsize=(2,2))

	if colors is None:
		colors=["blue"]*len(cov_arrays)


	
	relative_enrichments=[]
	for i in range(0,len(cov_arrays)):
		relative_enrichments.append(float(np.median(cov_arrays[i]))/(np.median(cov_ref_arrays[i])*CATCH_mass_fractions[i]))

	print(relative_enrichments)

	plt.bar(np.arange(len(relative_enrichments)), relative_enrichments, color=colors)
	plt.xticks(np.arange(len(labels)), labels)
	plt.ylabel("Relative enrichment")
	plt.show()

	if outpath is not None:
		plt.savefig(outpath, format="svg", bbox_inches="tight")
		plt.savefig(outpath.replace(".svg", ".pdf"), format="pdf", bbox_inches="tight")
	else:
		return fig
