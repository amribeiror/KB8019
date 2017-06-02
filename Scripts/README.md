## ORF Finder
* *find_orf_final.py*: Takes genome as arg[1] and find the possible open reading frames. User may need to edit script to specify whether they want to use transcription promoter filtering or translational promoter filtering or both. User will lso be asked whether the genome is a eukaryote or a prokaryote. The output is a *.orf_coords* file listing the start and stop positions of the ORFs and the corresponding reading frame.
* *extract_orflist_to_proteome.py*: Takes the genome as arg[1] and the *.orf_coords* as arg[2] for translation of the transcriptome to the proteome in a multifasta format.
* *orf_distribution_statistics*: This script should be run in the folder where the *.orf_coords* file are located. It uses globbing to aggregate all the files and construct a histogram to show the distribution of ORF length.

## Blastp Result Parser
* *xml_parser.py*: Takes blastp xml output as arg[1] and prints the statistics of only the best hits classified by e-value threshold of 10^-3.

## Genome and Proteome Statistics
* *_statistics.py*: Takes the whole genome file or a multifasta proteome file as arg[1] and calculate the statistics and write the them into a *.stats* file.
* *_matrix.py*: Takes the *.stats* file as arg[1] and create different distance matrices with each saved as a *.txt* file.
