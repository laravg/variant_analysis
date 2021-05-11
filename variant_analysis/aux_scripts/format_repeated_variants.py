import argparse
import os
from pathlib import Path
from itertools import zip_longest
import math
import pandas as pd
import numpy as np
from io import StringIO


def create_parser():
    parser = argparse.ArgumentParser(
        description='Gen 16S Oral Variants, all files will be generated in a folder in this directory')
    parser.add_argument('-p', '--repeated_variants_path', metavar='input', required=True,
                        help='path to the txt file with the repeated variants file')
    parser.add_argument('-f', '--fasta_path', metavar='input', required=True,
                        help='path to the txt file with the repeated variants file')
    parser.add_argument('-c', '--csv_path', metavar='input', required=True,
                        help='path to the txt file with the repeated variants file')
    return parser


def parse_arguments(parser):
    args = parser.parse_args()
    return args


def execute_parser():
    parser = create_parser()
    return parse_arguments(parser)


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

    # primero str
    # subsp
    # 2 palabras


class RepeatedVariant:
    header: str
    identifier: str
    sequence: str
    taxonomy_genomes: dict[str, str]
    new_header: str

    def __init__(self, header, sequence, seq_counter):
        self.header = header
        self.sequence = sequence
        self.identifier = f'RS{seq_counter:03d}'
        self.taxonomy_genomes = {}
        self.parse_header()

    def parse_header(self):
        splitted_header = self.header.split(', ')
        middle_position = math.floor(len(splitted_header)/2)
        middle = splitted_header[middle_position]
        first_ids = middle.split(' ')[-1]
        last_taxonomy = ' '.join(middle.split(' ')[:-1])
        splitted_header[middle_position] = last_taxonomy
        splitted_header.insert(middle_position + 1, first_ids)
        taxonomies = splitted_header[:len(splitted_header)//2]
        ids = splitted_header[len(splitted_header)//2:]
        for index, tax in enumerate(taxonomies):
            self.taxonomy_genomes[tax] = ids[index]
        self.new_header = f'>{self.identifier} {"/".join(self.taxonomy_genomes.values())}'

    def get_fasta_contents(self):
        return f'{self.new_header}\n{self.sequence}'

    def get_xlsx_contents(self):
        lines = []
        for index, ids in enumerate(self.taxonomy_genomes.values()):
            taxonomy = list(self.taxonomy_genomes.keys())[index]
            taxonomy_splitted = taxonomy.split(';')
            speciesNstrain = taxonomy_splitted[-1]
            species = ''
            strain = ''

            # if 'subsp.' in speciesNstrain:
            #     pos = speciesNstrain.find('subsp.')
            #     species = speciesNstrain[:pos]
            #     strain = speciesNstrain[pos:]
            # elif 'str.' in speciesNstrain:
            #     pos = speciesNstrain.find('str.')
            #     species = speciesNstrain[:pos]
            #     strain = speciesNstrain[pos:]
            # speciesNstrain_split = speciesNstrain.split(' ')
            # index = speciesNstrain_split.index('subsp.')
            # pos = index + 2
            # species_list = speciesNstrain_split[:pos]
            # species = ' '.join(species_list)
            # strain_list = speciesNstrain_split[pos:]
            # strain = ' '.join(strain_list)
            if 'sp.' in speciesNstrain and not 'subsp.' in speciesNstrain:
                if 'str.' in speciesNstrain:
                    pos = speciesNstrain.find('str.')
                    species = speciesNstrain[:pos]
                    strain = speciesNstrain[pos:]
                else:
                    species = speciesNstrain
                species = ' '.join(species.split(' ')[1:])

            else:
                speciesNstrain_split = speciesNstrain.split(' ')
                species_list = speciesNstrain_split[:2]
                species = ' '.join(species_list)
                strain_list = speciesNstrain_split[2:]
                strain = ' '.join(strain_list)
                species = ' '.join(species.split(' ')[1:])

            new_taxonomy = taxonomy_splitted[:-1] + \
                [species.strip()] + [strain.strip()]

            lines.append(
                f'{self.identifier};{ids};{";".join(new_taxonomy)}')
        return lines


args = execute_parser()

repeated_variants_path: str = args.repeated_variants_path
fasta_path: str = args.fasta_path
csv_path: str = args.csv_path

all_data = []
with open(repeated_variants_path, 'r', encoding="utf8") as repeated_variants_file:
    all_data = [line.strip() for line in repeated_variants_file.readlines()]

all_data_grouped = grouper(all_data, 2)

repeated_variants = []
seq_counter = 1
for header, seq in all_data_grouped:
    repeated_variants.append(RepeatedVariant(header, seq, seq_counter))
    seq_counter = seq_counter + 1

fasta_lines = [repeated_variant.get_fasta_contents()
               for repeated_variant in repeated_variants]
csv_header = ['Repeated_variant_id', 'Genome_ids', "Superkingdom", "Phylum", "Class", "Order",
              "Family", "Genus", "Species", "Strain"]
csv_contents = [repeated_variant.get_xlsx_contents()
                for repeated_variant in repeated_variants]
csv_lines = [';'.join(csv_header)] + \
    [item for sublist in csv_contents for item in sublist]

with open(fasta_path, 'w+') as fasta_file:
    fasta_file.write('\n'.join(fasta_lines))

csv_file = StringIO('\n'.join(csv_lines))

repeated_variants_df = pd.read_csv(csv_file, sep=';')

df_unique = repeated_variants_df.groupby('Repeated_variant_id')[[
    "Superkingdom", "Phylum", "Class", "Order",  "Family", "Genus", "Species", "Strain"]].nunique().add_prefix('num_').reset_index()

writer = pd.ExcelWriter(csv_path, engine='xlsxwriter')
workbook = writer.book

repeated_variants_df.to_excel(
    writer, sheet_name='repeated_variants', index=None, header=True)

worksheet = writer.sheets['repeated_variants']

# Get the dimensions of the dataframe.
(max_row, max_col) = repeated_variants_df.shape

# Create a list of column headers, to use in add_table().
column_settings = []
for header in repeated_variants_df.columns:
    column_settings.append({'header': header})

# Add the table.
worksheet.add_table(0, 0, max_row, max_col - 1,
                    {'columns': column_settings})

# Make the columns wider for clarity.
worksheet.set_column(0, max_col - 1, 20)


df_unique.to_excel(
    writer, sheet_name='counts', index=None, header=True)

worksheet = writer.sheets['counts']

# Get the dimensions of the dataframe.
(max_row, max_col) = df_unique.shape

# Create a list of column headers, to use in add_table().
column_settings = []
for header in df_unique.columns:
    column_settings.append({'header': header})

# Add the table.
worksheet.add_table(0, 0, max_row, max_col - 1,
                    {'columns': column_settings})

# Make the columns wider for clarity.
worksheet.set_column(0, max_col - 1, 20)

writer.save()
