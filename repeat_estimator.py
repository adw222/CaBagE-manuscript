"""
Updated 12-03-2020

Before starting:
minimap2 -Yax map-ont $REF $FASTQ > $FILE.sam
ml samtools
samtools view -Sb $FILE.sam > $FILE.bam
samtools sort -@ 8 -m 1G $FILE.bam -o $FILE.sorted.bam
samtools index $FILE.sorted.bam

Run script with following command:

repeat_estimator.py --bam --ref --locus --repeat_unit --alignment_buffer --id --allele-count

"""
import pysam
import argparse
from skbio.alignment import StripedSmithWaterman
from Bio import SeqIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import sklearn
from sklearn import mixture
from pyfaidx import Fasta

parser = argparse.ArgumentParser(description='find alignments that span a locus of interest')
parser.add_argument('--bam', dest='bam', help='The input bam file',
        required=True)
parser.add_argument('--ref', dest='ref', help='The indexed reference fasta file',
        required=True)
parser.add_argument('--locus', dest='locus', help='Postition of STR repeat in reference , e.g.: chr9:27573480-27573551', default = "chr9:27573485-27573546")
parser.add_argument('--repeat_unit', dest='repeat_unit', help='Repeated sequence of interest, e.g.: CCCCGG', default = "CCCCGG")
parser.add_argument('--alignment_buffer', dest='buf', default = "1000", help='buffer for alignment on either side of target, will vary by guide design')
parser.add_argument('--id', dest='id', default = "ND13803", help='Sample ID')
parser.add_argument('--allele_count', dest='ac', help='number of alleles present in indiviudal at locus of interest', default = "2")
args = parser.parse_args()


def get_spanning_read_ids(bamfilename, chrom, repeat_start, repeat_end, buf):
	"""
	Return a list of sequence ids having (possibly split) alignments
	that span the locus of interest.

	Assumption: reads that span should have a softclip alignment on the "left" of the repeat
	"""
	
	spanning_read_ids = {}
	#If file is BAM format:
	#bamfile = pysam.AlignmentFile(bamfilename, "rb")
	#If file is CRAM format:
	bamfile = pysam.AlignmentFile(bamfilename, "rc")
	# iterate over each alignment (aln) in the BAM file
	for aln in bamfile.fetch(chrom, repeat_start-buf, repeat_end+buf):
		# does the alignment "span" the repeat defined by start------end?
		# http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.infer_read_length
		if aln.reference_start <= repeat_start and aln.infer_read_length() > (repeat_end - aln.reference_start):
			spanning_read_ids[aln.query_name] = True
	#print(spanning_read_ids)
	print("NUMBER SPANNING READS:", len(spanning_read_ids))
	bamfile.close()
	return spanning_read_ids


def get_full_read_sequences(bamfilename, spanning_read_ids, chrom, repeat_start, repeat_end, buf):
	"""
	Return a mapping of query_name->full_dna_sequence
	for reads that span the locus.
	"""
	full_read_sequences = {}
	bamfile = pysam.AlignmentFile(bamfilename, "rb")
	for aln in bamfile.fetch(chrom, repeat_start-buf, repeat_end+buf):
		if aln.query_name in spanning_read_ids:
			if aln.query_name in full_read_sequences and aln.query_length >= len(full_read_sequences[aln.query_name]):
				full_read_sequences[aln.query_name] = aln.query_sequence
			elif aln.query_name not in full_read_sequences:
				full_read_sequences[aln.query_name] = aln.query_sequence

	bamfile.close()
	print("NUMBER OF SEQUENCES:", len(full_read_sequences))
	#for key, value in full_read_sequences.items() :
	#	print(key, value)
	return full_read_sequences

def estimate_repeat_number(target_sequences, reference, chrom, repeat_start, repeat_end):
	"""
	For each sequence that spans the locus, estimate the repeat 
	number from the distance between flanking unique sequence.
	"""
	ref = Fasta(reference)
	aln_lengths = []
	
	#query_sequence_A = "TTTTACTTAAAGAATTTCATCCACATCTTGTCAAGAGAGTTCAGTCTGATGGAAAGCACTGACTTCTATTTACAGAGCATTAGATGAGTGCTTTTATCATATTATGAGTAGGCATACAGAGCCTGGCAAAACAGTTAACTCTAAGTATGTACAGAAATGGTTGAACACAACGACAGTTTTAACACGTGTATTTGTAATTTCAAAAATTCATTTAGGTAATATTTACTTTTAAATATGTTGTATCAATTTAATAGTCTTAAGAGACAGCACTAGATATAAGCCGTACAGCTTCTTTAAAATATCCACTGTTTTTAATACAATGTAAGCAGTCAGTTTACAATGATCAAATATAGGAATGTAATCTGAATTGAAATGGTAATGACACTACTGCTGTCATAACTAACAACAGCAAACTGGAGGCCAACATAATGAATTAAGTTAACATACAACCATAAAATTATATTGCAAACATATTTTTCTTTCATTCTTTTAGGTTAAAAAGGTGGATAATCATAAAGGCAATATTACAACTCTAATATTTCATCATTAAACTGAAAATAAAAGTATTTCCTAAAACAGAACTGAACCCTGGAGCAAAATCTGATTGAATTATAGGGAAACTTTTACCACGTTGTGAAAATTGAACTATTATACTGCTAGTTACACTCTCACTCCTaacagaataagaaaaaaaaaatgggccgggcatggtgggtcacacctgttatcccagctctttggtaggccgaggcaggtggatcacctgaggtcaggagctcaagaccagcctggccaacatggtgaaaccccacctctactaaaaatacaaaaaattagccgggtgtggtggtggacacctgtaatcccagctactcgggaggctgaggcaggagaatctcttgaacccggaggtggcagaggttgcaatgagctgagatggcgccactgcactccagcctgggcgacagagagagactctgcctcaagaaaaaaacaaacaaacaaacaaacaaaaagaataagaaagaaaatgaaGGACAAAGATCATACTGAATTGCTTAGTTTTAAATCCTACCAAAAGAAATAGCCTGGGAAATGAAATGTCACAGAGAAGTATAATCAGGAGAGCTGTACAATTATTTTACTAATACTTGAAGTCATCGTCTTTGGTGAGAAAAATCCATACATGCAAATGCAGCTGAAAAAAATCAGCTCAAAACCAATAGTTGTTTATGTACCTATCTTACGTACATGTAGTGCTGTCTACTCCAGAGAGTTACCAAACATTAGCCAGTCTTTTGAGGGAAGCCAAGATTCAAATTGAGTGAGACGGTGGCTTGCTCACAGGGTTCATGAGAGGTTTCCCAATACACTTTCTGGAAATAATCCcatacatgcagacatgattacattaattaacatctgctaaaactgttagtagagtgctaagtttgaggttttgctttttctttaaacgtctgttaaaaaatcaaccatctcttccctgattggtatttagaaaggtggttggtccactgctattgtAGTGAAAATTCTACAATCATAAAGCCCTCACttcttgttttttagagacagggtctcgttttgtcatccaggctggaatgcactggcaggatcatagctctcggtaacttcaaactcttgggctcaaatgaccctcctgcctcagcctcccaagtagctaggactacaggtgcacatcaccacgcccggctaagtttttaattttttgtagagacagggtctacgttgcccaggttgagcttgaactcctggcttcaagtgatcctcttgcctccgcctcccaaagctctggcattacaggtgtaagccaccTTCTCCAACCTGGCTCTCAATACTTGTAACCATGCTGTTTATTTTCTCCCAGCCCAAAGAGAAGCAGGATCCTAAACCGTCCACTTTCCACAACAGGAGCTGCCCAGGACCACTTCAAGGACAGTGAACTGTTTACAGTACCAGAAAGTTCACAACACTTTCTCAATCTTCAACATCAGGGAAGACTGGAAGGTGAAGTTCATATCACTATCTGGCCATTTCTCACAGTTCCAAGTTTCTCAGACAATAGGTAGGCTAACCTAGTCCTCCTGGGAACTATCTAATTAACGTAGAATAGAACCCGAGGGCAGACTTGAAAAACAGAAGTCCTCCTTGGTTTACTTTGTTTCTCTGAAAGCAAATTGTGGAGTGCCAACATAGCCAAACAAAATATTTTATCAACTTCATAAGGTGCTTGTAATTTTTTCCTGGAGCAGGTAAATGCTGGCTTAGTGAACAATCTGGAATGTGGTAATTACTCTCGTTCTTGTTTCAGATGTACTATCAGCATGTAGCAGTTTCCAACTGATTCAGGGTTTTCCTAAAGTGGCAGGCCTTGGCAGAGGTGGTGACAACAATGCCCGTGTCAAATGACACCGTATTTCAAGTATTCTGACTCCAGGTTATTAATATCCcctatatgatagtcttgtttctgtgatattcacagattatgttaaaagtttcccaaagtctgagaaaaatcatatcttaacagtatcttttttttttttgatcctttgtacaaaagtagaagtaatgccagacagattacgtacccttgttgtgaacaactggtgcatggcaactgtttgaatagaaatttaccaactgccacaaccaggcaactactctcccagagcctaacaatctcgatTGCATCTGAAAGGGCCACCCCTCCTGGGAAAGTGCAGGACCTCCCTCCTGTTTCTGAATACAAAGCCTGGTGGTGTTCAACGCGGCCAGATAGACCCAATGAGCACACGGACATGTAATCTGTGCACTTCTTTAGACAACTGATTACCATCAGTCAAGTGATGCCCAAGTCACAATAGTCACTTCCTTTAAGCAAGTCTGTGTCATCTCGGAGCTGTGAAGCAACCAGGTCATGTCCCACAGAATGGGGAGCACACCGACTTGCATTGCTGCCCTCATATGCAAGTCATCACCACTCTCTAGAAGCTTGGGCTGAAATTGTGCAGGCGTCTCCACACCCCCATCTCATCCCGCATGATCTCCTCGCCGGCAGGGACCGTCTCGGGTTCCTAGCGAACCCCGACTTGGTCCGCAGAAGCCGCGCGCCGCCCACCCTCCGGCCTTCCCCCAGGCGAGGCCTCTCAGTACCCGAGGCTCCCTTTTCTCGAGCCCGCAGCGGCAGCGCTCCCAGCGGGTCCCCGGGAAGGAGACAGCTCGGGTACTGAGGGCGGGAAAGCAAGGAAGAGGCCAGATCCCCATCCCTTGTCCCTGcgccgccgccgccgccgccgccgccgGGAAGCCCGGGGCCCGGATGCAGGCAATTCCACCAGTCGCTAGAGGCGAAAGCCCGACACCCAGCTTCGGTCAGAGAAATGAGAGGGAAAGTAAAAATGCGTCGAGCTCTGAGGAGAGCCCCCGCTTCTACCCGCGCCTCTTCCCGGCAGCCGAACCCCAAACAGCCACCCGCCAGGATGCCGCCTCCTCACTCACCCACTCGCCACCGCCTGCGCCTCCGCCGCCGCGGGCGCAGGCACCGCAACCGCAgccccgccccgggcccgcccccgggcccgccccgaccacg"
	#query_sequence_B = "ccccTAGCGCGCGACTCCTGAGTTCCAGAGCTTGCTACAGGCTGCGGTTGTTTCCCTCCTTGTTTTCTTCTGGTTAATCTTTATCAGGTCTTTTCTTGTTCACCCTCAGCGAGTACTGTGAGAGCAAGTAGTGGGGAGAGAGGGTGGGAAAAACAAAAACACACACCTCCTAAACCCACACCTGCTCTTGCTAGACCCCGCCCCCAAAAGAGAAGCAACCGGGCAGCAGGGACGGCTGACACACCAAGCGTCATCTTTTACGTGGGCGGAACTTGTCGCTGTTTGACGCACCTCTCTTTCCTAGCGGGACACCGTAGGTTACGTCTGTCTGTTTTCTATGTGCGATGACGTTTTCTCACGAGGCTAGCGAAATGGGGCGGGGCAACTTGTCCTGTTCTTTTATCTTAAGACCCGCTCTGGAGGAGCGTTGGCGCAATAGCGTGTGCGAACCTTAATAGGGGAGGCTGCTGGATCTGGAGAAAGTGAAGACGATTTCGTGGTTTTGAATGGTTTTGTTTgtgcttggtaggcagtgggcgctcaacacataattggtggatgaaATTTTGTTTTTACCGTAAGACACTGTTAAGTGCATTCAAAACTCCACTGCAAACCCTGGTAGGGGACAGCTCCGGCACTGCGGGCGGGAATCCCACGGTCCCCTGCAAAGTCATCGCAATTTTGCCTTTACATGTAAGAATTCTCTCAAGCATGATTTTCACACTGGGGAATGTCATTTTTGCTAGTTGCAATATGTGGATGAGTTGTTTTTTTTTAACTTTTGAAAAACGTACCATTCTGTTTGATGTGTAAAAAACACAAAGATTTTTGAAACCTTGCGTCTTTTGGTCTGCAGGTGTATAGATTCCACTTACTACAGATGAGTAGCATTTACACCACTCAGATGTGTAAAAAAACAAAGGTTTTTTAAACTGTGTGCCTTTTGATCTGCAAGTGTGAGATGGCACTTACTACAGTGAGTAGCATTTAATCTTTTTCATCACTAAAAATCACACAGAACGTTTTAATCATTCACCGAGGAAGAAAGGGAGGAATAAATACACAAAATGGCTCTCAACGTCTACACCTTCTGCAGAAACAGACCCTTTTCCTACTGTTCTATGCTTtgtgaaagttgatcatacaaattgggtcattctttttatacccaactaaaatagtgggggtagggggtagaaaagcacttaggacaaatgacactgctcccacagtgtaattctctccaagtccagctgctgcaactgcccgttgtgacctgagaccagttttatctaatagttgctaaaatgacctgctgcagctctaattttatctaccaccatcactcaccagttgaaactcaccagctcctcagatccttaatagtgccaatgaattttctcaaagagcactatgtaacatttctcttttttaacaaaacctcccccttttctttgttgtgtggatataccgaagaccatctgatctacatgtatgccctaattgcaattctttcttcccaaataaatcacttaatttagagattcatctctgtatttttattttgactgacaGCTTATAACAAGTAGCTAGCATTTACCAAGTTTCTACACTGAGTTGTACTtcacttatacgtggaattaaaaaacaactgaatttatagaaacagagtagacccttggttggggggcttggggggaaagaaaattgtagggtagggtacaaagttgcagttacgtctaatacatctagagatttaatgtacaacatgaggactagcgttaataattgtgttagtccattcttacactgctataaagaaataactgaaactgggtaatttataaagaaaagtttaatggctcacagttctgcaggctgtacaagaagcatggctggatcagcttctgggcaggccatagggaacttaaaatcatgatggaaggcatagggagaccccagacttcacatggcaggaactgggggaagagagaaatgggaggtgctacatacgtttaaacaactagatcttgtcagaactcactatatagtaccaagaggggactgtacaaaaccattagaagccaccccataatccactcacctcccaccaggcccaacctccaacactggggattacagttgaacatgagatttgggtggggacagagatccaaaccatgttattccaactctggcccctcccaaatctaatgtccttctcatattgcaaaatactgtcgtgccttaccaacagttccccaaagtcttaactcgatccagcattcattcaaaagtccaaagtcccaagtctcacctgagacgaagctagtcccttctacctatgaacctgtaaaatcaaaaacaaggtaattgcttcaaagatacaatgggggtataggcattgggcagatactgccattccgaaagggagaaatctgccaaaagaaagaggctatagggccccattgcaagtctgaaagccagccgggcagtcattaaatgttaaagctctgaaataatctcctttgactcacacccagggaacactgatgcaatgagtgggctcccaaaaccttgggcagaaccacccctgtggttttccagggttcatctcccacagctgctctcatgggctagcattgagtgcttgcagcttttccaggctgcagggtgcaagttgttggtggatctaccattctggggtctggaggacggtggctgtcttgtcatagctctgctaggcagtgccccaggggactctctgtgggggctgcaaccccacatttcttctccttgcttccctagtagatgttctccatgaggattccaccccagtaacaggcttctgtctggacatccaggctttttcatacatcctctaaaatctaggcagagcttcttaagcctcaactcttgcattatgtgcgcccgccggcttcacagcttatggaagccaccaaggcttatgcctggcaccctgtgaagcagcagcctgaactgtattcttactggtgaaagttatctgagttaccagctgcaaatccatgtgggtctgcagcaacctcaattcttgcctcctcagaagaaagaatttgaccaagaggcataaggcagaaaaagagactgcgacaagtttcagagcaggagtaaaagtttattaaaaagctttagaacaggaatgaaaggaaagtacatttggaagaggcccaagtgggcaccttggaggtcaagtgccctgtttgaccttgaacctaggatcttatacactggcctacttctgacatcttgtgcccctttcccttggtccttccctaagggtgagcttgccgcatgcatggtgccctgcttgcacttggaaggtgagcgtgtgcagtgtgtttactggagttgtatacatgcttacctgaggctttcttcccttttccggtggaatgcccccaaaggtcatacttcaccattttgcctcttaatgtgcatgttaagcccactctctcagttcctgagatcttattggaagcgcccagttaccaatttcaggtgtttctatctattgagaagttgcctctccctggtgctggctgcaaccaattactatttt"
	query_sequence_A = ref[chrom][(repeat_start-1000):(repeat_start)]
	query_sequence_B = ref[chrom][(repeat_end):(repeat_end+1000)]
	#print(query_sequence_A)
	queryA = StripedSmithWaterman(str(query_sequence_A), gap_open_penalty=2, gap_extend_penalty=1)
	queryB = StripedSmithWaterman(str(query_sequence_B), gap_open_penalty=2, gap_extend_penalty=1)
	for qname, target_sequence in target_sequences.items():
		alignmentA = queryA(target_sequence)
		alignmentB = queryB(target_sequence)

		if int(alignmentB.target_begin) > int(alignmentA.target_end_optimal):
			a_len = [qname,alignmentA.optimal_alignment_score, alignmentB.optimal_alignment_score, (alignmentB.target_begin - alignmentA.target_end_optimal+1), alignmentA.target_begin, alignmentA.target_end_optimal, alignmentB.target_begin, alignmentB.target_end_optimal]
			aln_lengths.append(a_len)

	print("NUMBER ALIGNMENT LENGTHS:", len(aln_lengths))
	return aln_lengths

def rep_stats(bamfilename, chrom, repeat_start, repeat_end, buf, alignment_lengths, id, repeat_unit, ploidy):
	#Takes list of lists
	data = pd.DataFrame(alignment_lengths, columns = ['readID','alnAscore', 'alnBscore', 'alnLEN', 'alnAs', 'alnAe', 'alnBs', 'alnBe'])

	bamfile = pysam.AlignmentFile(bamfilename, "rb")
	strands = {}
	strand = []
	for aln in bamfile.fetch(chrom, repeat_start-buf, repeat_end+buf):
		strands[aln.query_name] = aln.is_reverse
	for read in data['readID']:
		strand.append(strands.get(read))

	data['strand'] = strand

	ct = pd.crosstab(data['strand'], data['alnLEN']>100)
	#print(ct)
	c2 = scipy.stats.chi2_contingency(ct, correction=True, lambda_=None)
	#print(c2)
	#MANN WHITNEY TEST
	#MW = scipy.stats.mstats.mannwhitneyu(data['strand'], data['alnLEN'], use_continuity=True)
	#print(MW)
	#tt = scipy.stats.ttest_ind(data['strand'], data['alnLEN'], axis=0, equal_var=True, nan_policy='propagate')
	#print(tt)

	fig = plt.figure()
	kwargs = dict(histtype='stepfilled', alpha=0.3, bins=40)

	unit_length = len(repeat_unit)
	plot_title = id + " Repeat Copy Number Distribution"

	plt.hist(data['alnLEN']/unit_length, bins=40) 

	plt.xlabel('Repeat Copy Number')
	plt.ylabel('Read Count')
	plt.title(plot_title) 
	fig.savefig("%s.pdf" % id, transparent=True)

	# Fitting using a GMM with a single component
	clf = mixture.GaussianMixture(n_components=int(ploidy)) ###DEFAULT PLOIDY 2, CAN CHANGE AFTER LOOKING AT HISTOGRAM IN CASE OF MOSAICISM
	v = np.array(data['alnLEN'])
	emd = v.reshape(-1, 1)
	clf.fit(emd)
	print("Repeat copy number (means, weights): ", clf.means_/unit_length, clf.weights_)

	#plt.show()
	#print(data)



def main():
	"""
	This function drives the script's order of execution
	"""
	(chrom, start, end) = args.locus.replace('-', ':').split(':')
	start = int(start)
	end = int(end)
	buf = int(args.buf)
	#print(args.bam)
	#print(start)

	# Step 1, get a dictionary of read ids of spanning sequences
	spanning_read_ids = get_spanning_read_ids(args.bam, chrom, start, end, buf)
	
	# Step 2, get a dictionary of the full sequence of these ids
	full_read_sequences = get_full_read_sequences(args.bam, spanning_read_ids, chrom, start, end, buf)

#	fastq_count = estimate_repeat_number_fastq(args.fastq)

	# Step 3.
	alignments = estimate_repeat_number(full_read_sequences, args.ref, chrom, start, end)

#	rep_stats(alignments, fastq_count, args.id)
	rep_stats(args.bam, chrom, start, end, buf, alignments, args.id, args.repeat_unit, args.ac)
	

if __name__ == "__main__":
    main()


