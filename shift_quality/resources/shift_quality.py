
import sys
import scipy.stats
import numpy as np
import random
import argparse
from itertools import izip

# This function converts the sum of all quality values into an integer value
def calculate_quality_score(quality_string, input_file, item_name):

    error_probability = 0
    for x in quality_string:
        try:
            phred = ord(x)
            if phred < 33 or phred > 75:
                raise ValueError()
            error_probability += 10**-((phred-33)/10.0)
        except ValueError:
            # If a quality value outside of Illumina1.8 format is detected, exit with an informatice message
            sys.stderr.write("ERROR: Encountered a non-supported quality value - %s\n" % x)
            sys.stderr.write("Only Illumina1.8 format quality values are supported\n")
            sys.stderr.write("Found in the %s input file in read %s - %s\n" % (input_file, item_name, quality_string) )
            sys.exit(1)

    return error_probability


# This function measures the properites of the reads in this chunk of the file and writes those reads which are consistent with the shifted
# distribution of qualities specified
def extract_read(read_info, mate_info, error_probability, downsample_fraction, original_distribution, shifted_distribution):


    # Calculating the exponential probabilities very far from the mean causes floating point errors in numpy. Capping the maxium errors
    # prevents calculating extremely small probabilities in typical quality distributions. In practice, for any shift a read with this value
    # will have an inclusion probability of either 0.0 or 1.0
    if error_probability > 50.0:
        error_probability = 50.0

    # The probability to include should be a function of how many excess reads we have, and how often we would see a read of this quality
    # in the desired distribution as compated to the starting one
    inclusion_probability = downsample_fraction * shifted_distribution.pdf(error_probability) / original_distribution.pdf(error_probability)

    if inclusion_probability > random.random():
        reads_output_file.write("".join(read_info))
        mates_output_file.write("".join(mate_info))
        return True

    return False

parser = argparse.ArgumentParser(description='Extract reads to shift a starting distribution of qualities')

parser.add_argument('--reads', action="store", dest="reads", required=True)
parser.add_argument('--mates', action="store", dest="mates", required=True)
parser.add_argument('--output-reads', action="store", dest="output_reads", required=True)
parser.add_argument('--output-mates', action="store", dest="output_mates", required=True)
parser.add_argument('--downsample-fraction', action="store", dest="downsample_fraction", type=float, default=1.0)
parser.add_argument('--starting-mean', action="store", dest="starting_mean", type=float, required=True)
parser.add_argument('--starting-stdev', action="store", dest="starting_stdev", type=float, required=True)
parser.add_argument('--stdev-shift', action="store", dest="stdev_shift", type=float, required=True)

arguments = parser.parse_args()

downsample_fraction = arguments.downsample_fraction
standard_deviation_shift = arguments.stdev_shift
mean_quality = arguments.starting_mean
stdev_quality = arguments.starting_stdev

reads_output_file = open(arguments.output_reads, 'w')
mates_output_file = open(arguments.output_mates, 'w')

# Line numbers will be used to determine which of the 4 lines per fastq read we are reading
line_number = 0

# Record how many pairs we have looked at
total_pairs = 0
total_extracted = 0

# Record the qualities included in batches
qualities = []

# Precalculate the distributions for the current quality distribution and for the desired one. This will be used to determine 
# the probability to extract a read to shift the distribution
original_distribution = scipy.stats.norm(mean_quality, stdev_quality)
shifted_distribution = scipy.stats.norm(mean_quality + standard_deviation_shift * stdev_quality, stdev_quality)

# Loop over the reads and the mates together
with open(arguments.reads) as reads_file, open(arguments.mates) as mates_file: 
    for line, mate_line in izip(reads_file, mates_file):

        # The first line is the read name, currently this doesn't check to make sure the read names match
        if line_number % 4 == 0:
            read_name = line
            mate_name = mate_line

        # The second line is the nucleotide sequence
        if line_number % 4 == 1:
            read_nt = line
            mate_nt = mate_line

        # The third line is the spacer: +
        # The fourth line contains the read qualities, which is what we will use to determine the overall quality of this read pair
        if line_number % 4 == 3:
            read_info = [read_name, read_nt, "+\n", line]
            mate_info = [mate_name, mate_nt, "+\n", mate_line]

            # Make a unified score for this read based on the sum of the phred values of all qualities in the read and mate (assumes Illumina1.8 format)
            error_probability =  calculate_quality_score(line.strip(), arguments.reads, read_name)
            error_probability += calculate_quality_score(mate_line.strip(), arguments.mates, mate_name)

            # Do a biased random selection based on the desired final distribution relative to the current one
            # Write the read to the output and return whether it was selected
            is_extracted = extract_read( read_info, mate_info, error_probability, downsample_fraction, original_distribution, shifted_distribution )

            total_pairs += 1

            # Collect and output statistics on the number and distribution of the extracted reads every 100,000 extracted reads
            # Note if you change the buffer size that a large buffer will become the determinant of memory use by the script
            if is_extracted:
                total_extracted += 1
                qualities.append(error_probability)
                if total_extracted % 100000 == 0:
                    extracted_mean = np.mean(qualities)
                    extracted_stdev = np.std(qualities)
                    sys.stderr.write("Processed %d reads, Extracted %d reads, Mean quality this 100K: %.2f +- %.2f\n" % (total_pairs, total_extracted, extracted_mean, extracted_stdev))
                    qualities = []

        line_number += 1

reads_file.close()
mates_file.close()
reads_output_file.close()
mates_output_file.close()
