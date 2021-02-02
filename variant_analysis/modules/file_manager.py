from pathlib import Path
from Bio import SeqIO
import os
import classes.Genome as Genome
from modules import downloads


def create_dir(path):
    if not os.path.isdir(path):
        Path(path).mkdir(parents=True, exist_ok=True)


def read_input_dir(path):
    # TODO: check if genomes repeated in file, if so, skip and warn
    genomes = []
    genome_ids = []
    for filename in os.listdir(path):
        if os.path.isfile(f'{path}/{filename}') and (filename.endswith(".fasta") or filename.endswith(".fa")):
            fasta_sequences = SeqIO.parse(
                open(f'{path}/{filename}'), 'fasta')
            for fasta in fasta_sequences:
                header, identifier, sequence = fasta.description, fasta.id, str(
                    fasta.seq)
                genome = Genome(header, identifier, filename, sequence)
                if identifier in genome_ids:
                    raise Exception('Repeated identifier, ignoring')
                else:
                    genome_ids.append(identifier)
                    genomes.append(genome)

    return genomes


def read_input_file(path):
    # TODO: check if ids repeated in file, if so, skip and warn
    genome_ids = []
    if path.endswith(".txt") or path.endswith(".csv"):
        with path.open('r') as file:
            ids = file.read().splitlines()
            genome_ids = [id.strip() for id in ids if id]

    return genome_ids
