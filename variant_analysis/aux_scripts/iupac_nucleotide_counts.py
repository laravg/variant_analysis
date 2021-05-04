import argparse
import os
from pathlib import Path
from typing import Optional
import json
from Bio import Entrez
from Bio import SeqIO
import xml.etree.ElementTree as ET
from classes.Genome import Genome
import re


SLEEP_BETWEEN_TRIES = 5
API_KEY = '8eee4fbae5b295870dd6c2669300f1334008'  # Must use one
EMAIL = 'laramaria.vazquez.gonzalez@rai.usc.es'
TOOL = 'gen16sov'

Entrez.email = EMAIL
Entrez.api_key = API_KEY
Entrez.tool = TOOL
Entrez.sleep_between_tries = SLEEP_BETWEEN_TRIES


class Genome:

    header: str
    identifier: str
    sequence: str
    filename: str

    def __init__(self, header, identifier, sequence, filename):
        self.header = header
        self.identifier = identifier
        self.sequence = sequence
        self.filename = filename

    def get_fasta_contents(self):
        return f'{self.header}\n{self.sequence}'


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
    parser.add_argument('-n', '--nucleotides', metavar='nucleotides', nargs='+', required=True,
                        help='nucleotides iupac')
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


def read_input_file(path: str) -> tuple[list[str], Optional[str]]:
    genome_ids = []
    unique_genome_ids = []
    msg = None
    with open(path, 'r') as file:
        lines = file.read().splitlines()
        genome_ids = [id.strip()
                      for line in lines if line for id in line.split(',')]
        unique_genome_ids = list(set(genome_ids))
        duplicate_genome_ids = list(set(
            [x for x in genome_ids if genome_ids.count(x) > 1]))
        if duplicate_genome_ids:
            msg = f'Ignored repeated genome ids: {duplicate_genome_ids}'
    if not unique_genome_ids:
        raise Exception('No genomes found in file')
    return (unique_genome_ids, msg)


def read_input_dir(path: str) -> tuple[list[Genome], Optional[str]]:
    genomes = []
    genome_ids = []
    msg = None
    for filename in os.listdir(path):
        if os.path.isfile(f'{path}/{filename}') and (filename.endswith(".fasta") or filename.endswith(".fa")):
            fasta_sequences = SeqIO.parse(
                open(f'{path}/{filename}'), 'fasta')
            for fasta in fasta_sequences:
                header, identifier, sequence = fasta.description, fasta.id, str(
                    fasta.seq)
                genome = Genome(f'>{header}', identifier,
                                sequence, None)
                if identifier not in genome_ids:
                    genome_ids.append(identifier)
                    genomes.append(genome)

    unique_genome_ids = list(set(genome_ids))
    if len(unique_genome_ids) != len(genome_ids):
        msg = f'Ignored repeated genome ids: {[identifier for identifier in genome_ids if identifier not in unique_genome_ids]}'

    if not genomes:
        raise Exception('No genomes found in directory')
    return (genomes, msg)


def get_genomes(genome_ids: list[str]) -> list[Genome]:
    handle = Entrez.efetch(
        db="nucleotide", id=genome_ids, rettype="fasta", retmode='xml')
    genomes = __parse_entrez_list(Entrez.read(handle))
    handle.close()

    return genomes


def __parse_entrez_list(entrez_list) -> list[Genome]:
    genomes = []
    for elem in entrez_list:
        genome_id = str(elem['TSeq_accver'])
        description = str(elem['TSeq_defline'])
        header = f'>{genome_id} {description}'
        sequence = str(elem['TSeq_sequence'])
        genome = Genome(header, genome_id, sequence, None)
        genomes.append(genome)
        with open(f'use_case/downloaded_genomes/{genome_id}.fa', 'w+') as genome:
            genome.write(header + '\n' + sequence)

    return genomes


def save_iupac_nucleotide_counts(genomes: list[Genome], nucleotides, path: str):

    iupac_codes = ['N', 'B', 'V', 'M', 'Y', 'S', 'K', 'R', 'W', 'H', 'D']
    iupac_codes_header = [code + '>10' for code in iupac_codes]

    iupac_nucleotide_counts_header = ['Genome']
    # iupac_nucleotide_counts_header.extend(nucleotides)
    # iupac_nucleotide_counts_header.extend(iupac_codes_header)
    iupac_nucleotide_counts_header.append('Any')
    lines = [';'.join(iupac_nucleotide_counts_header)]

    with open(path, 'w+') as iupac_nucleotide_counts:
        iupac_nucleotide_counts.write(lines[0] + '\n')
        for genome in genomes:
            line = [genome.identifier]
            # for nuc in nucleotides:
            #     line.append(str(genome.sequence.count(nuc)))
            # for code in iupac_codes:
            #     pattern = code + '{11,}'
            #     line.append(str(len(re.findall(pattern, genome.sequence))))
            pattern = '[NBVMYSKRWHD][NBVMYSKRWHD]'
            distinct = 0
            for find in re.findall(pattern, genome.sequence):
                if find[0] != find[1]:
                    distinct = distinct + 1
            line.append(str(distinct))
            iupac_nucleotide_counts.write(';'.join(line) + '\n')


args = execute_parser()

input_path: str = args.input
output_path: str = args.output
nucleotides: list[str] = args.nucleotides


genomes = []

if os.path.isfile(input_path):
    print('reading')
    genome_ids, msg = read_input_file(input_path)
    print('downloading')
    genomes = get_genomes(genome_ids)
elif os.path.isdir(input_path):
    print('reading')
    genomes, msg = read_input_dir(input_path)

print('counting')
save_iupac_nucleotide_counts(genomes, nucleotides, 'countsDEFbact2.csv')
