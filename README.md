# HTT Trinucleotide Repeat Characterization
# [![DOI](https://zenodo.org/badge/838543003.svg)](https://doi.org/10.5281/zenodo.14861905)

This set of scripts is a workflow for characterizing the CAG and CCG repeats in exon1 of HTT.

This workflow was created for and used in data analysis for the GEM-HD paper: "Genetic modifiers of somatic expansion and clinical phenotypes in Huntington’s disease reveal shared and tissue-specific effects" (preprint: https://www.biorxiv.org/content/10.1101/2024.06.10.597797v1; doi: https://doi.org/10.1101/2024.06.10.597797).

## Workflow Overview

### **htt_trinucleotide_characterization.py**

Will read through flattened fastq files and analyze read(pairs) for the CAG and CCG repeat sequences in HTT. Will begin characterizing from the first instance of 3 CAG repeats until a sequence of CAGCTTCCT is encountered (downstream of the CCG repeat) or the end of the read is reached.

Parameters

[REQUIRED]

sample: "Sample name."

miseq_run: "Name of this miseq run."

fastq: "File with reads, flattened."

output: "Output prefix."

quality: "Quality System used during sequence.", choices=['p33', 'p64']

[OPTIONAL]

--fastq2: "If supplied, will perform paired end analysis. Paired-files must be sorted in the same order and contain only valid pairs. Program will terminate if an incorrect pair is encountered."

--graph":: "Indicate whether to produce graphical results"

--test_orient: "Whether to test the orientation of each read sequence"

--prolines: "Will perform Prolines-Only Characterization from nononucleotides backwards."

--cag_only: "Will seek to only characterize CAG sequence from 3 CAGs to a stop of 3 CCGs"

--qstats: "Will write out quality scores statistics for the sequence. Currently only supporting proline analysis."



### **miseq_allelic_extractor.py**

Takes the CAG_DISTRIBUTION output file and RESULTS output file of **htt_trinucleotide_characterization.py** and will analyze the distribution of CAG sizes and abundant sequences to determine the CAG genotype and corresponding phased sequence.

sample: "Sample name."

cag_dist: "CAG Distribution of the MiSeq analyzed sample."

res_f: "Result file with structure frequencies."



### **miseq_instability.py**

Takes the CAG_DISTRIBUTION output file and the TOP_ALLELES output file and will analyze the instability of the expanded allele.

[REQUIRED]

cag: "Miseq-analyzed CAG distribution."

alleles: "Identified CAG alleles for this sample."

[OPTIONAL]

--peak_bias": "Peak-bias for when to begin summing peak-proportions. Default = 0 (begin at main peak)

--cutoff_override: "Override 35 CAG cutoff."

--a1: "Allele 1 manual override."

--a2: "Allele 2 manual override."

## EXAMPLE

See the example folder for prepared input data and expected output data for your verification and validation purposes of running these scripts.



[![DOI](https://zenodo.org/badge/838543003.svg)](https://doi.org/10.5281/zenodo.14861905)

