#!/usr/bin/env python

import argparse as ap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.axes as maxes

plt.switch_backend("agg")

__author__ = 'kc199'

NONONUCLEOTIDE = "CAGCTTCCT"


class Sequence:
    def __init__(self, fastq_record, qual, orientation_test=False):
        self.id = fastq_record[0]
        seq = fastq_record[1]
        if orientation_test:
            # Compare ratio of CAG to CTG present in sequence. If CAG/CTG < 1, reverse complement
            cag_count = seq.count("CAG")
            ctg_count = seq.count("CTG")
            try:
                if cag_count / ctg_count < 1:
                    self.seq = dna_reverse_complement(seq)
                else:
                    self.seq = seq
            except ZeroDivisionError:
                self.seq = seq
        else:
            self.seq = seq

        self.qual = fastq_record[3]
        self.qual_system = qual
        self.qscores = self.qual_symbols_to_scores(self.qual)

    def findcag(self, strict=False):
        # Find three CAG to avoid matching a CAG that may be upstream of the repeat region and any missenses
        # if strict, demand upstream 6 bases of CAG to also be present
        if strict:
            pos = self.seq.find("AAGTCCTTCCAGCAGCAG")
            if pos == -1:
                return pos
            else:
                return pos + 9
        else:
            return self.seq.find("CAGCAGCAG")

    def find_nononucleotide(self):
        # Find the nononucleotide (CAGCTTCCT)
        return self.seq.find(NONONUCLEOTIDE)

    def create_flat_fastq(self):
        return [self.id, self.seq, '+', self.qual]

    def qual_symbols_to_scores(self, quality_symbols):
        try:
            if self.qual_system == "p33":
                return np.array([ord(q) - 33 for q in quality_symbols])
            elif self.qual_system == "p64":
                return np.array([ord(q) - 64 for q in quality_symbols])
        except:
            return np.array([0] * len(quality_symbols))


def dna_reverse_complement(s):
    dna_complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'N': 'N'
    }

    rna_complement = {
        'A': 'U',
        'U': 'A',
        'C': 'G',
        'G': 'C',
        'N': 'N'
    }
    return "".join([dna_complement[n] for n in s])[::-1]





def process_read(seq: Sequence):
    cag_start = seq.findcag()
    if cag_start == -1:
        unprocessed.append(seq)
        return 0, 0, 0
    else:
        triplet = 1
        triplet_type = 'CAG'
        triplet_count = 0
        trip_seq = []

        while True:
            iter_start = cag_start + (triplet - 1) * 3
            iter_stop = cag_start + (triplet * 3)

            if iter_stop > len(seq.seq):
                trip_seq.append([triplet_type, triplet_count])
                quality = seq.qscores[cag_start: iter_stop - 3]
                break
            elif seq.seq[iter_start: iter_stop + 6] == NONONUCLEOTIDE:
                trip_seq.append([triplet_type, triplet_count])
                trip_seq.append([NONONUCLEOTIDE, '1'])
                quality = seq.qscores[cag_start: iter_stop + 6]
                break
            else:
                triplet += 1
                trip = seq.seq[iter_start: iter_stop]
                if trip != triplet_type:
                    trip_seq.append([triplet_type, triplet_count])
                    triplet_type = trip
                    triplet_count = 1
                else:
                    triplet_count += 1

    if trip_seq[0][0] == "CAG":
        cag_count = trip_seq[0][1]
    else:
        cag_count = 0

    return trip_seq, cag_count, quality


def process_reads_cag(seq: Sequence):
    cag_start = seq.findcag(strict=True)
    if cag_start == -1:
        unprocessed.append(seq)
        return 0, 0, 0
    else:
        triplet = 1
        triplet_type = 'CAG'
        triplet_count = 0
        trip_seq = []

        while True:
            iter_start = cag_start + (triplet - 1) * 3
            iter_stop = cag_start + (triplet * 3)

            if iter_stop > len(seq.seq):
                trip_seq.append([triplet_type, triplet_count])
                quality = seq.qscores[cag_start: iter_stop - 3]
                break
            elif seq.seq[iter_start: iter_stop + 6] == 'CCGCCGCCG':
                trip_seq.append([triplet_type, triplet_count])
                trip_seq.append(['CCGCCGCCG', '1'])
                quality = seq.qscores[cag_start: iter_stop + 6]
                break
            else:
                triplet += 1
                trip = seq.seq[iter_start: iter_stop]
                if trip != triplet_type:
                    trip_seq.append([triplet_type, triplet_count])
                    triplet_type = trip
                    triplet_count = 1
                else:
                    triplet_count += 1

    if trip_seq[0][0] == "CAG":
        cag_count = trip_seq[0][1]
    else:
        cag_count = 0

    return trip_seq, cag_count, quality


def process_reads_proline(seq: Sequence):
    nono_start = seq.find_nononucleotide()
    if nono_start == -1:
        unprocessed.append(seq)
        return 0, None, None
    else:
        triplet = 1
        triplet_type = NONONUCLEOTIDE
        triplet_count = 1
        trip_seq = []

        while True:
            iter_start = nono_start - (triplet - 1) * 3
            iter_stop = nono_start - (triplet * 3)

            if iter_stop <= 0:
                trip_seq.append([triplet_type, triplet_count])
                quality = seq.qscores[iter_stop + 3: nono_start + 9]
                break
            elif seq.seq[iter_stop - 6: iter_start] == "CAGCAGCAG":
                trip_seq.append([triplet_type, triplet_count])
                trip_seq.append(["CAGCAGCAG", 1])
                quality = seq.qscores[iter_stop - 6: nono_start + 9]
                break
            else:
                triplet += 1
                trip = seq.seq[iter_stop: iter_start]
                if trip != triplet_type:
                    trip_seq.append([triplet_type, triplet_count])
                    triplet_type = trip
                    triplet_count = 1
                else:
                    triplet_count += 1

    return list(reversed(trip_seq)), None, quality


def aggregate_structure(struct, qscores=None, proline_quality_file=None):
    # Hash the structure and put in dictionary. Increment already existing structure by one or make new
    struct_index = hash(str(struct))
    try:
        aggregate[struct_index]["count"] += 1
        if qscores is not None:
            aggregate[struct_index]["quality"].append(qscores)
    except KeyError:
        aggregate[struct_index] = {"structure": struct, "count": 1, "quality": [qscores]}
    finally:
        if qscores is not None and proline_quality_file is not None:
            record = [structure_serializer(struct), str(struct_index), " ".join([str(q) for q in qscores])]
            proline_quality_file.write("\t".join(record) + '\n')


def print_outputs(agg, cag, unproc, graph_output, mismatch=None, quality=None, paired=True):
    if paired:
        filename_prefix = "%s_%s" % (args.output, "paired_consensus")
    else:
        filename_prefix = args.output

    result = []
    for k, v in agg.items():
        structure = v["structure"]
        structure = [f'{trip}{qty}' for trip, qty in structure]
        structure = " ".join(structure)
        count = v["count"]
        if quality:
            median_quality = np.median(np.array(v["quality"]), axis=0)
            result.append([args.sample, structure, median_quality, count])
        else:
            result.append([args.sample, structure, count])

    if not args.prolines:
        with open("%s.cag_distribution" % filename_prefix, 'w') as cag_dist:
            for c in cag:
                cag_dist.write("%s\n" % c)

    if quality:
        cols = ['sample', 'structure', 'median_quality', 'count']
    else:
        cols = ['sample', 'structure', 'count']

    result_df = pd.DataFrame(result, columns=cols)
    result_df.sort_values(by="count", ascending=False, inplace=True)
    if not result_df.empty:
        result_df['freq'] = result_df['count'] / result_df['count'].sum()

    result_df["miseq_run"] = args.miseq_run
    
    result_df.to_csv("%s.results.txt" % filename_prefix, sep='\t', index=False)

    with open("%s.unprocessed_reads" % filename_prefix, 'w') as un:
        for s in unproc:
            un.write("\t".join(s.create_flat_fastq()) + '\n')

    if mismatch:
        with open("%s.mismatched_pairs" % filename_prefix, 'w') as mm:
            for m in mismatch:
                mm.write("\t".join(m) + '\n')

    if graph_output and not result_df.empty and not args.prolines:
        viz(result_df, cags, filename_prefix)

    # Print analysis statistics
    with open("%s.statistics" % filename_prefix, 'w') as stats:
        header = "\t".join(['sample', 'total', 'processed', 'unprocessed', 'mismatched'])
        stats.write(header + '\n')
        proc = len(cag)
        nonproc = len(unproc)
        subtotal = proc + nonproc
        if mismatch:
            mism = len(mismatch)
            info = "\t".join([str(x) for x in [args.output, subtotal + mism, proc, nonproc, mism]])
        else:
            info = "\t".join([str(x) for x in [args.output, subtotal, proc, nonproc]])

        stats.write(info + '\n')


def print_qscore_output(agg):
    with open("%s.quality_score_stats" % args.output, 'w') as qstats:
        header = ["structure", "structure_count", "structure_sequence", "min_quality_sequence", "median_quality_sequence",
                  "max_quality_sequence"]
        qstats.write("\t".join(header) + '\n')
        results = []
        for k, v in agg.items():
            median_quality = [str(x) for x in np.median(np.array(v["quality"]), axis=0)]
            max_quality = [str(x) for x in np.max(np.array(v["quality"]), axis=0)]
            min_quality = [str(x) for x in np.min(np.array(v["quality"]), axis=0)]

            raw_struct = v['structure']
            structure = structure_serializer(v['structure'])

            structure_sequence = "".join([x[0] * int(x[-1]) for x in raw_struct])

            results.append([structure, str(v['count']), structure_sequence, " ".join(min_quality), " ".join(median_quality), " ".join(max_quality)])

        results.sort(key=lambda x: int(x[1]), reverse=True)
        for r in results:
            qstats.write("\t".join(r) + '\n')


def viz(agg, cag, filename_prefix):
    fig, (ax, axtable) = plt.subplots(ncols=1, nrows=2, figsize=(9.6, 5.4))
    # Plot CAG histogram distribution
    bins = np.arange(0, max(cag) + 10)
    assert isinstance(ax, maxes._axes.Axes)
    ax.hist(cag, bins=bins)
    ax.set_ylabel("Number of Sequences")
    ax.set_xlabel("CAG Length")
    ax.set_title("%s - CAG Distribution" % args.output)

    # Make table of top characterizations. Limit to frequency > 0.01
    table_data = agg[agg.freq > 0.01]
    if len(table_data) > 10:
        table_data = table_data[:10]
    table_data['structure'] = table_data.apply(func=lambda x: structure_serializer(x['structure']), axis=1)
    pd.options.display.max_colwidth = 100
    table_string = table_data.to_string(header=False, index=False)

    axtable.text(0, 0, table_string, fontsize=10, horizontalalignment='left')
    axtable.axis('off')

    fig.savefig("%s.graphs.svg" % filename_prefix)


def structure_serializer(x):
    a = ["".join([str(j) for j in i]) for i in x]
    return ", ".join(a)


parser = ap.ArgumentParser(description="This script will read in a fastq file of reads (or a pair of fastqs for PE) "
                                       "and analyze the structure of the CAG region in HTT and produce the top most "
                                       "frequent structures and a CAG distribution.")

parser.add_argument("sample", help="Sample name.")
parser.add_argument("miseq_run", help="Name of this miseq run.")
parser.add_argument("fastq", help="File with reads, flattened.")
parser.add_argument("output", help="Output prefix.")
parser.add_argument("quality", help="Quality System used during sequence.", choices=['p33', 'p64'])
parser.add_argument("--fastq2", help="If supplied, will perform paired end analysis. Paired-files must be sorted "
                                     "in the same order and contain only valid pairs. Program will terminate "
                                     "if an incorrect pair is encountered.")
parser.add_argument("--graph", action="store_true", help="Indicate whether to produce graphical results")
parser.add_argument("--test_orient", action="store_true", help="Whether to test the orientation of each read sequence")
parser.add_argument("--prolines", action="store_true", help="Will perform Prolines-Only Characterization from"
                                                            " nononucleotides backwards.")
parser.add_argument("--cag_only", action="store_true", help="Will seek to only characterize CAG sequence from 3 CAGs to"
                                                            " a stop of 3 CCGs")
parser.add_argument("--qstats", action="store_true", help="Will write out quality scores statistics for the sequence. "
                                                          "Currently only supporting proline analysis.")
args = parser.parse_args()
if args.prolines:
    args.output = "%s_%s" % (args.output, "Proline")

if args.cag_only:
    args.output = "%s_%s" % (args.output, "cagonly")

# Process reads
unprocessed = []
aggregate = {}
cags = []
if args.fastq2:
    mismatch_structure = []
    with open(args.fastq, 'r') as fq, open(args.fastq2, 'r') as fq2:
        for line1, line2 in zip(fq, fq2):
            if line1.strip() == 'empty' or line2.strip() == 'empty':
                break

            record1, record2 = line1.strip().split("\t"), line2.strip().split('\t')
            if record1[0].split(" ")[0] != record2[0].split(" ")[0]:
                raise SystemExit("Nonmatching FastQ read pair encountered.")
            else:
                fastq1, fastq2 = Sequence(record1, args.quality, args.test_orient), Sequence(record2, args.quality,
                                                                                             args.test_orient)
                if args.prolines:
                    structure1, cag1, qscores1 = process_reads_proline(fastq1)
                    structure2, cag2, qscores2 = process_reads_proline(fastq2)

                elif args.cag_only:
                    structure1, cag1, qscores1 = process_reads_cag(fastq1)
                    structure2, cag2, qscores2 = process_reads_cag(fastq2)
                else:
                    structure1, cag1, qscores1 = process_read(fastq1)
                    structure2, cag2, qscores2 = process_read(fastq2)

                if structure1 != structure2:
                    mismatch_structure += [fastq1.create_flat_fastq() + fastq2.create_flat_fastq()]
                elif structure1 == structure2 == 0:
                    continue
                else:
                    if cag1 is not None:
                        cags.append(cag1)

                    aggregate_structure(structure1)

    print_outputs(aggregate, cags, unprocessed, args.graph, mismatch=mismatch_structure)

else:
    with open(args.fastq, 'r') as fq:
        if args.prolines and args.qstats:
            qsf = open("%s.quality_sequences" % args.output, 'w')
        # Read each fastq record
        for line in fq:
            record = line.strip().split('\t')
            fastq = Sequence(record, args.quality, args.test_orient)
            # Analyze the structure, return 0 if no CAG start site found
            if args.prolines:
                structure, cag, qscores = process_reads_proline(fastq)
            elif args.cag_only:
                structure, cag, qscores = process_reads_cag(fastq)
            else:
                structure, cag, qscores = process_read(fastq)
            if structure != 0:
                if cag is not None:
                    cags.append(cag)
                if args.prolines and args.qstats:
                        aggregate_structure(structure, qscores=qscores, proline_quality_file=qsf)
                else:
                    aggregate_structure(structure, qscores=qscores)
            else:
                continue

        try:
            qsf.close()
        except NameError:
            pass

    if args.prolines:
        if args.qstats:
            print_qscore_output(aggregate)
        print_outputs(aggregate, cags, unprocessed, args.graph, paired=False)
    else:
        print_outputs(aggregate, cags, unprocessed, args.graph, paired=False)
