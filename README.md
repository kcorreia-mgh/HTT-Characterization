# HTT-Characterization

This set of scripts is a workflow for characterizing the CAG and CCG repeats in exon1 of HTT.

**htt_trinucleotide_characterization.py**
Will read through flattened fastq files and analyze read(pairs) for the CAG and CCG repeat sequences in HTT. Will begin characterizing from the first instance of 3 CAG repeats until a sequence of CAGCTTCCT is encountered (downstream of the CCG repeat) or the end of the read is reached.
