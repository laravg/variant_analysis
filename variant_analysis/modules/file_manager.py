from pathlib import Path
from typing import Optional
from Bio import SeqIO
import os
from classes.Genome import Genome
# import csv


def create_dir(path: str):
    if not os.path.isdir(path):
        Path(path).mkdir(parents=True, exist_ok=True)


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
                                None, sequence, filename)
                if identifier not in genome_ids:
                    genome_ids.append(identifier)
                    genomes.append(genome)

    unique_genome_ids = list(set(genome_ids))
    if len(unique_genome_ids) != len(genome_ids):
        msg = f'Ignored {len(genome_ids)-len(unique_genome_ids)} repeated genomes'

    if not genomes:
        raise Exception('No genomes found in directory')
    return (genomes, msg)


def read_input_file(path: str) -> tuple[list[str], Optional[str]]:
    genome_ids = []
    unique_genome_ids = []
    msg = None
    with open(path, 'r') as file:
        lines = file.read().splitlines()
        genome_ids = [id.strip()
                      for line in lines if line for id in line.split(',')]
        unique_genome_ids = list(set(genome_ids))
        if len(unique_genome_ids) != len(genome_ids):
            msg = f'Ignored {len(genome_ids)-len(unique_genome_ids)} repeated genomes'
        # elif '.csv' == path.suffix:
        #     sniffer = csv.Sniffer()
        #     first_lines = file.read(1024)
        #     has_header = sniffer.has_header(first_lines)
        #     dialect = sniffer.sniff(first_lines)
        #     file.seek(0)
        #     csvreader = csv.reader(file, dialect)
        #     if has_header:
        #         next(csvreader, None)
        #     lines = list(csvreader)
        #     print(dialect.delimiter)
        #     genome_ids = [id.strip()
        #                   for line in lines if line for id in line.split(dialect.delimiter)]
    if not unique_genome_ids:
        raise Exception('No genomes found in file')
    return (unique_genome_ids, msg)
