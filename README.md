# Hsum
Script to document runs of homozygosity in genome-wide microarray geontype data.

This script will read a genome-wide genotype data file and output a list of runs of homozygosity.  
It will read through all of the text files in the indicated directory, and generate one output file per input.
The input file must be a tab delimited text file with columns: marker name, dbSNP rs number, chromosome, physical position, genetic position.  Genotype calls must be AA, AB, BB, or NO.
Command line options:
 --help, print this informaiton and quits.
 -t, number of threads to use for processing.
 -p, % of non-homozygous SNPs allowed in homozygous block, default 5.
 -h, max number of heterozygous SNPs allowed in a homozygous block, default 1000 (not used)
 -c, max number of consecutive heterozygous SNPs allowed, default 2.
 -x, number of SNPs at the end of a block to monitor for excess heterozygous snps, default 10
 -y, number of heterozygous SNPs allowed in the monitored SNPs at the end of a block, default 3
 -1, first criteria to use for sorting output, can be 'markers', 'bp', or 'cM', default cm.
 -2, second criteria to use for sorting output, must be one of the options above, and cannot be the same given for -1, default bp.
 --cmcut, cutoff used for smallest cM block reported, default 0.5.
 --bpcut, cutoff used for smallest bp block reported, default 1.
 --markcut, cutoff used for smallest number of markers in a block, default 25
 --dir, locaton of the genotype files, default “\.”
