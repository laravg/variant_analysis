import argparse
import os
from pathlib import Path


def execute_parser():
    parser = create_parser()
    return parse_arguments(parser)


def create_parser():
    parser = argparse.ArgumentParser(
        description='Gen 16S Oral Variants, all files will be generated in a folder in this directory')
    parser.add_argument('-i', '--input', metavar='input', required=True, type=valid_input,
                        help='path to the txt file with the Genbank IDs or directory with FASTA files to analyze')
    parser.add_argument('-o', '--output', metavar='output', required=True, type=valid_output,
                        help='path to the directory where the results will be stored')
    parser.add_argument('-f', '--format', metavar='format', nargs='+', required=True, choices=['mothur', 'qiime', 'dada2'],
                        help='output formats: mothur, qiime, dada2 (separated by spaces)')
    parser.add_argument('-td', '--training-database', metavar='training_database', choices=['HOMD', 'GREENGENES', 'ARCHAEA', 'ESCAPA'], default='HOMD',
                        help='HOMD for Human Oral Microbiome Database 16S rRNA (default), GREENGENES for Greengenes 16S rRNA gene database')
    parser.add_argument('-tl', '--taxonomy-level', metavar='taxonomy_level', nargs='+', choices=["superkingdom", "phylum", "class", "order", "family", "genus", "species"],
                        default='species', help='Taxonomy level')  # TODO: change help
    parser.add_argument('-s', '--substitute_iupac', type=valid_output,
                        help='path to folder to store genomes after substitute IUPAC characteres')
    parser.add_argument('-p', '--primers', metavar='primers', type=valid_input,
                        help='path to the csv file with the primer sequences (forward and reverse) to trim the found genes')
    return parser


def parse_arguments(parser):
    args = parser.parse_args()
    return args


def valid_input(param):
    if not os.path.isfile(param) and not os.path.isdir(param):
        raise argparse.ArgumentTypeError('File/directory not accessible')
    if os.path.isfile(param) and not Path(param).suffix == '.txt' and not Path(param).suffix == '.csv':
        raise argparse.ArgumentTypeError('File must have .txt extension')
    return param


def valid_output(param):
    if os.path.isdir(param) and os.listdir(param):
        raise argparse.ArgumentTypeError(f'Directory {param} not empty')
    return param
