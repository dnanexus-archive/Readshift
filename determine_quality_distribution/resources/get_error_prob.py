import sys
import numpy as np
import argparse
from itertools import izip
import gzip

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

parser = argparse.ArgumentParser(description='Extract reads to shift a starting distribution of qualities')

parser.add_argument('--reads', action="store", dest="reads", required=True)
parser.add_argument('--mates', action="store", dest="mates", required=True)

# The quality value is kept in memory until the end of the program. For large numbers of reads this can use significant amounts
# of memory. This argument causes the program to exit after hitting this limit.
parser.add_argument('--maximum-reads', help="Compute statistics on at most this many reads", action="store", dest="maximum_reads", required=False)

parser.add_argument('--output-quality-score-file', action="store", dest="output_quality_score_file", required=True)

arguments = parser.parse_args()

# Record the number of lines iterated over for
line_number = 0
bases = 0

# Loop over the reads and the mates together
output_quality_score = gzip.open(arguments.output_quality_score_file, 'wb')
with gzip.open(arguments.reads, 'rb') as reads_file, gzip.open(arguments.mates, 'rb') as mates_file:
    for line, mate_line in izip(reads_file, mates_file):

        # The first line is the read name, currently this doesn't check to make sure the read names match
        if line_number % 4 == 0:
            read_name = line
            mate_name = mate_line

        # The fourth line contains the read qualities, which is what we will use to determine the overall quality of this read pair
        if line_number % 4 == 3:

            # Make a unified score for this read based on the sum of the phred values of all qualities in the read and mate (assumes Illumina1.8 format)
            error_probability =  calculate_quality_score(line.strip(), arguments.reads, read_name)
            error_probability += calculate_quality_score(mate_line.strip(), arguments.mates, mate_name)

            # Record the number of bases for coverage calculations
            bases = len(line.strip()) + len(mate_line.strip())

            # Append the value seen to an array. The size of this will determine the ultimate memory requirements
            # it is better for performance to pass only a subset of the FASTQ file to this.
            output_quality_score.write('{}\t{}\n'.format(error_probability, bases))

        line_number += 1

        # Check if we have exceeded the total number of reads to process
        if arguments.maximum_reads is not None:
            if line_number > int(arguments.maximum_reads) * 4:
                break

output_quality_score.close()
