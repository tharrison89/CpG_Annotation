"""
Entry point for the annotator workflow
"""

import argparse
import gzip
import os
from Bio import SeqIO
import sys

from annotator.kmeans_gc_content_annotator import KmeansGcContentAnnotator

OPENING_FUNCTIONS = {True: gzip.open, False: open}
EXPECTED_GZ_EXTENSION = ".gz"

FASTA_FILE_TYPE = "fasta"
FASTQ_FILE_TYPE = "fastq"
SEQ_FILE_TYPES = {FASTQ_FILE_TYPE, FASTA_FILE_TYPE}

class UnsupportedFileException(Exception):
    """
    Raised when an unsupported file type is used as input
    """
    pass

def parse_seq_file(input_seq_file, seq_file_type):
    """
    Iterate records of a sequencing file

    :param str input_seq_file: Path to input sequence file
    :param str seq_file_type: The sequence filetype
    :yields str, str: Sequence record id and sequence
    """
    if seq_file_type not in SEQ_FILE_TYPES:
        msg = "{} is not a supported input file type".format(seq_file_type)
        raise UnsupportedFileException(msg)

    is_gzed = input_seq_file.endswith(EXPECTED_GZ_EXTENSION)

    # a different opening function will be called based on whether the file has a gz extension
    with OPENING_FUNCTIONS[is_gzed](input_seq_file) as input_fh:
        for record in SeqIO.parse(input_fh, seq_file_type):
            yield record.id, str(record.seq).upper()

def run_annotator(sequence_file, file_type, output_file_handle,
                  kmer_length=KmeansGcContentAnnotator.DEFAULT_KMER_LENGTH,
                  min_feature_size=KmeansGcContentAnnotator.DEFAULT_MIN_FEATURE_SIZE,
                  alpha=KmeansGcContentAnnotator.DEFAULT_ALPHA):
    """
    Run the annotator

    :param str sequence_file: Path to the sequence file
    :param str file_type: input file type see SEQ_FILE_TYPES for available options
    :param file output_file_handle: Output file handle
    :param int kmer_length: Desired Kmer length
    :param int min_feature_size: Minimum feature size
    :param float alpha: Alpha for significance
    """

    annotator = KmeansGcContentAnnotator(min_feature_size, kmer_length)

    for record_id, sequence in parse_seq_file(sequence_file, file_type):
        output_file_handle.write(">{}{}".format(record_id, os.linesep))
        intervals = annotator.annotate_sequence(sequence)
        for interval in annotator.iterate_significant_hits(alpha):
            output_file_handle.write("{}{}".format(interval, os.linesep))

def parse_arguments():
    """
    Parse command line args
    :return Namespace args: namespace from the arg parser
    """
    parser = argparse.ArgumentParser(description='Annotate a sequence by GC content.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', "--seq_file", type=str, required=True,
                        help="input sequence file to annotate")

    parser.add_argument('-k', "--kmer_length", type=int, required=False,
                        default=KmeansGcContentAnnotator.DEFAULT_KMER_LENGTH,
                        help="Sliding window length")

    parser.add_argument('-m', "--min_interval_length", type=int, required=False,
                        default=KmeansGcContentAnnotator.DEFAULT_MIN_FEATURE_SIZE,
                        help="Minimum length of an interval")

    parser.add_argument('-t', "--input_file_type", type=str, required=False,
                        default=FASTA_FILE_TYPE, choices=SEQ_FILE_TYPES,
                        help="input file format")

    parser.add_argument('-a', "--alpha", type=float, required=False,
                        default=KmeansGcContentAnnotator.DEFAULT_ALPHA,
                        help="Alpha for significance")

    parser.add_argument('-o', "--outfile", type=argparse.FileType("w"), required=False,
                        default=sys.stdout,
                        help="output file. Default is to write to stdout")

    args = parser.parse_args()
    return args

def main():
    """
    Main method
    """

    args = parse_arguments()

    run_annotator(args.seq_file, args.input_file_type, args.outfile,
                  kmer_length=args.kmer_length,
                  min_feature_size=args.min_interval_length,
                  alpha=args.alpha)

if __name__ == "__main__":
    main()
