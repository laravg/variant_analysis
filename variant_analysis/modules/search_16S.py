"""
search_16S.py - a python implementation of the 16S identification method described by
Robert Edgar. See preprint: https://www.biorxiv.org/content/biorxiv/early/2017/04/04/124131.full.pdf

For the kmer counting, I've used khmer's [https://github.com/dib-lab/khmer] CountGraph, which works
well for small values of k (like the default, 13), but obviously gets memory expensive with larger k.

Usage:
    search_16S.py train -i REF_SEQS -d REF_DB [-v] [-k KMER_SIZE]
    search_16S.py find -q SEQS -d REF_DB [-v] [-m MAX_MISMATCH -o OUTPUT --min_len MIN_LEN --max_len MAX_LEN]
    search_16S.py train find -i REF_SEQs -q SEQS -d REF_DB [-v] [-k KMER_SIZE -m MAX_MISMATCH -o OUTPUT --min_len MIN_LEN --max_len MAX_LEN]

Options:
    -i REF_SEQS           FASTA of known high quality 16S sequences. R.Edgar's preprint uses Greengeens v13.5 97% OTUs
    -q SEQS               FASTA of contigs to be searched.
    -d REF_DB             Filename for saving/loading the kmer countgraph. This is not a human readable file
    -k KMER_SIZE          Size of sliding windows used to tokenize reference sequence [default: 13]
    -m MAX_MISMATCH       Maximum number of mismatches allowed in the start or end motif for 16S search [default: 4]
    -o OUTPUT             File to write identified 16S sequences and fragments to. Defaults to stdout
    --min_len MIN_LEN     Minimum length of reported sequences [default: 1000]
    --max_len MAX_LEN     Maximum length of reported sequences [default: 2500]
    -v, --verbose         Print progress stats to stderr
    -h, --help            Print this help text


"""
import io
from io import StringIO
import logging
import sys
import os

from Bio import SeqIO
import khmer
import numpy as np
import regex  # Note: this is regex, not re, because we want fuzzy matching
import skbio
from skbio import DNA
from docopt import docopt

histogram_window = 1000
initial_motif = "GNTTGATCNTGNC"
final_motif = "AGTCNNAACAAGGTANCNNTA"


def train(ref_seqs, kmer_size, ref_db):
    nkmers = 4 ** kmer_size
    tablesize = nkmers + 10
    kmer_counts = khmer.Countgraph(kmer_size, tablesize, 1)
    seqs_consumed, kmers_consumed = kmer_counts.consume_seqfile(ref_seqs)
    # logging.info('Parsed {} reference sequences which had {} kmers'.format(
    #     seqs_consumed, kmers_consumed))
    kmer_counts.save(ref_db)


def find(sample_seqs, ref_db, outfile_handle, max_n_mismatch_in_motif, min_len, max_len):
    kmer_counts = khmer.load_countgraph(ref_db)
    kmer_size = kmer_counts.ksize()
    chance_of_kmer_appearing = kmer_counts.n_unique_kmers()/4.00**kmer_size
    # logging.info("Proportion of reference kmers to all possible kmers: {}".format(
    #     chance_of_kmer_appearing))

    if '>' in sample_seqs:
        sample_seqs = StringIO(sample_seqs)
    seqs_to_assess = skbio.io.read(sample_seqs, format='fasta')
    motif = regex.Regex(
        '(?b)(?:{}){{s<={}}}(?:[ACGT]{{{},{}}})(?:{}){{s<={}}}'.format(skbio.DNA(initial_motif).to_regex().pattern,
                                                                       max_n_mismatch_in_motif,
                                                                       min_len,
                                                                       max_len,
                                                                       skbio.DNA(
                                                                           final_motif).to_regex().pattern,
                                                                       max_n_mismatch_in_motif))

    nseqs_processed = 1
    for seq in seqs_to_assess:
        length_seqs = []
        genes = []
        output_header = []
        output_gene_16s = []

        for i in range(2):
            strand = '+'
            if i == 1:
                strand = '-'
                seq = DNA(seq, positional_metadata={
                          'quality': range(len(seq))})
                seq = seq.reverse_complement()

            found_kmers = np.array(kmer_counts.get_kmer_counts(str(seq))) > 0
            windowed_sums = found_kmers.cumsum()
            windowed_sums[histogram_window:] = windowed_sums[histogram_window:] - \
                windowed_sums[:-histogram_window]
            # Finding stretches of high density, using code from
            # https://stackoverflow.com/questions/38161606/find-the-start-position-of-the-longest-sequence-of-1s
            candidate_segments = np.where(np.diff(np.hstack(([False],
                                                             windowed_sums > (
                                                                 histogram_window/2),
                                                             [False]))))[0].reshape(-1, 2)
            candidate_segments[:, 1] += kmer_size - 1
            candidate_segments[:, 0] -= int((histogram_window/2) + 200)
            candidate_segments[:, 1] += int((histogram_window/2) + 200)
            candidate_segments[candidate_segments[:, 0] < 0, 0] = 0
            candidate_segments[candidate_segments[:, 1]
                               >= len(seq), 1] = len(seq) - 1

            for pair in candidate_segments:
                if (pair[1] - pair[0]) < min_len:
                    continue
                subseq = seq[np.arange(*pair)]
                # Search for whole motif
                whole_match = regex.search(motif, str(subseq))

                if whole_match:
                    output_header.append('>' + seq.metadata['id'] + ";" +
                                         np.array2string(pair[0] + whole_match.span(), separator=",") + strand + ';\n')
                    output_gene_16s.append(whole_match.captures()[0])
                    length_seqs.append(len(whole_match.captures()[0]))

            # if nseqs_processed % 100 == 0:
            #     logging.info(
            #         'Processed {} test sequences'.format(nseqs_processed))
            nseqs_processed += 1

        bad_gene = None

        if length_seqs:
            max_length = max(length_seqs)
            min_length = min(length_seqs)
            diff_length = max_length - min_length

            if diff_length > 12:
                bad_gene = length_seqs.index(max_length)

        for index in range(len(output_gene_16s)):

            if bad_gene is None or (bad_gene is not None and index is not bad_gene):
                outfile_handle.write(output_header[index])
                outfile_handle.write(output_gene_16s[index] + '\n')


def main(args):
    if args["train"]:
        train(args["-i"], int(args["-k"]), args["-d"])
    if args["find"]:
        if args["-o"] is None:
            args["-o"] = sys.stdout
        else:
            args["-o"] = open(args["-o"], "w+")

        find(args["-q"], args["-d"], args["-o"], int(args["-m"]),
             int(args["--min_len"]), int(args["--max_len"]))

        if args["-o"] is not sys.stdout:
            args["-o"].close()


if __name__ == '__main__':
    args = docopt(__doc__)
    if args["--verbose"]:
        logging.basicConfig(level=logging.INFO)
        logging.info(args)
    main(args)
