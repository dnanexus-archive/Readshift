# Readshift

Readshift is a method to generate real NGS datasets of a desired quality profile through biased downsampling
from a dataset of larger coverage. 

This can be used on high coverage, well characterized genomes to create datasets that emulate low quality sequencing
runs. It can also be used to emulate the quality profile of a given sequencing project on a benchmark set to quantify
the False Positive and False Negative rates likely to be seen in the analysis of the dataset.

Readshift operates on pairs of reads. It calculates the expectation value for the total number of errors in a read+mate
by converting the Phred-encoded base qualtilies to probabilities and taking their sum. For each read+mate pair, 
Readshift calcualtes at what probability the pair should be included in order to reach a specified shifted distribution,
based on a specified starting distribution.

Readshift does not simulate data, nor does it alter any bases or quality values. It simply downsamples selectively.

## Evaluation Data

We have applied Readshift to 200X Hiseq2500, 200X HiseqX, and 320X Novaseq data, shifting each by 0.5, 1.0, 1.5, and 2.0
standard deviations in quality. This data is available from the DNAnexus platform. You will need an account to access it
but there is no charge to creating an account, accessing, or downloading this data.

The FASTQ and VCFs from running these through GATK3, GATK4, Sentieon, DeepVariant, Strelka2, and Freebayes pipelines can
be accessed at [this DNAnexus project titled Readshift](https://platform.dnanexus.com/projects/F9K5zYQ0ZQBb5jJf6X2zZPv1/data/)

An investigation on this data is detailed in our [Readshift blog post](https://blog.dnanexus.com/2018-01-16-evaluating-the-performance-of-ngs-pipelines-on-noisy-wgs-data/)

## Dependencies

Readshift requires the numpy and scipy Python libraries. 

## Inputs

Readshift operates on paired FASTQ files. It requires both the forward and reverse reads in uncompressed format. To apply
Readshift to compressed files, we suggest the use of names pipes. An example is present in the **shift_quality/shift_quality.sh** 
script which the DNAnexus app uses. Readshift expects that the FASTQ encoding format is Illumina1.8+ and will raise an exception
if it detects quality values outside of this range.

Small test inputs to demonstrate the scripts are found in the **test_data** directory. 

## Running Readshift

This repository contains two scripts required to apply Readshift:

**get_quality_distribution.py** - is used to profile fastq files. It outputs average fastq quality, standard deviation in
quality, and coverage. 

**usage: python get_quality_distribution.py --reads test_data/reads.R1.fastq --mates test_data/reads.R2.fastq**

**shift_quality.py** - is used to generate new FASTQ files with a new quality distribution. You will need to specify how to
model this shift - the mean and standard deviation of the starting FASTQ and the desired standard deviation shift

**python shift_quality.py --reads test_data/reads.R1.fastq --mates test_data/reads.R2.fastq --output-reads reads.shifted.R1.fastq --output-mates reads.shifted.R2.fastq --starting-mean 6.7 --starting-stdev 3.2 --stdev-shift 0.5 --downsample-fraction 0.5**

This will apply a shift of 0.5 standard deviations based on a starting normal distribution of 6.7 errors per read+mate
and a standard deviation of 3.2 errors. It will target a downsample fraction of 0.5 (downsample fraction of 0.0 loses all reads, 1.0
tries not to apply any downsample aside from the shift).

You may see a warning: **RuntimeWarning: invalid value encountered in double_scalars** - if this occurs, it means that numpy
is calculating exceedingly small probabilities of seeing a given read error profile based on the distributions. If this occurs,
you may want to cap the maximum error number in the extract_read function to a smaller maximum value.

## DNAnexus apps

The folders **determine_quality_distribution** and **shift_quality** contain the code for platform apps on DNAnexus. These show
how I apply the two scripts in an efficient way in a shell scripts. The apps reference an associated asset which contains the
pre-installed numpy and scipy dependencies. If you are curious how apps on DNAnexus work, these are examples of a few, simple
ways of building them.

Assuming you have installed the DNAnexus SDK and have logged into the platform, you may build them with 
**dx build determine_quality_distribution** and **dx build shift_quality**

## Notes about shifts and distributions

Readshift models the starting and target distribution of read+mate expected errors as normally distributed. In practice,
it seems that the real distribution of errors is a left-leaning curve. Because the use of the normal distribution only
serves to model the transformation, the final distribution will retain the shape of the original and will not be normally
distributed.

However, the output coverage of the shift will be less than the downsample fraction, because you must discard reads that are
a better fit to the old distribution than to the new one in order to apply a shift. The magnitude of possible shift that can
be achieved is a function of total coverage.

## Runtime

For a full 50X WGS run, **get_quality_distribution.py** and **shift_quality.py** each take about a day to run. I have not made
specific effort to optimize for speed. Most of the sequencing events for my initial investigation were from many, small FASTQ files
that were fast to run multiple files in parallel.

## Areas for improvement

I will happily consider pull requests and suggestions to improve Readshift. A number of areas for improvement come to mind:

**Performance** - Simply using multiprocessing should probably allow a speedup of 4x-6x on machines with enough cores. 

**Distributions** - Using a different model from the normal distribution might allow for more interesting transformations. 
In particular, it could be interesting to identify shift machineries which emphasize including many slightly worse reads 
as opposed to those with more signficant issues (or vice-versa).

**Error Modeling** - The error profile of reads is a function of their sequence context. It would be good to be able to 
have a correction factor which helped adjust reads that come from difficult regions from those which suffer from problems
during sequencing. This likely won't make a large impact in positive shifts, however, I think that negative shifts will be
limited in utility until this is corrected. With a high-quality, high-coverage dataset there are not enough bad reads to
apply a significant negative distribution shift (that is, trying to make the data cleaner) without starting to lose good
reads that come from regions you want to capture.


