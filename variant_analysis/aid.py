import pandas as pd
import numpy as np
import argparse

categories_order = {"superkingdom": 0, "phylum": 1, "class": 2, "order": 3,
                    "family": 4, "genus": 5, "species": 6}


def create_parser():
    parser = argparse.ArgumentParser(
        description='Aid')
    parser.add_argument('-d', '--directory', metavar='input', required=True,
                        help='path')
    return parser


def parse_arguments(parser):
    args = parser.parse_args()
    return args


def execute_parser():
    parser = create_parser()
    return parse_arguments(parser)


def concat_genomes(genomes):
    genomes = list(filter(lambda id: id != '_', genomes))
    if genomes:
        return '-'.join(genomes)
    else:
        return '__'


def separate_species(speciesNstrain):
    species = ''
    strain = ''

    if 'sp.' in speciesNstrain and not 'subsp.' in speciesNstrain:

        if 'str.' in speciesNstrain:
            pos = speciesNstrain.find('str.')
            species = speciesNstrain[:pos]
            strain = speciesNstrain[pos:]
        else:
            species = speciesNstrain

        species = ' '.join(species.split(' ')[1:])
        #species.split(' ')

    else:
        speciesNstrain_split = speciesNstrain.split(' ')
        species_list = speciesNstrain_split[:2]
        species = ' '.join(species_list)
        strain_list = speciesNstrain_split[2:]
        strain = ' '.join(strain_list)
        species = ' '.join(species.split(' ')[1:])

    return species.strip()


def separate_strain(speciesNstrain):
    species = ''
    strain = ''

    if 'sp.' in speciesNstrain and not 'subsp.' in speciesNstrain:
        if 'str.' in speciesNstrain:
            pos = speciesNstrain.find('str.')
            species = speciesNstrain[:pos]
            strain = speciesNstrain[pos:]
        species = speciesNstrain

    else:
        speciesNstrain_split = speciesNstrain.split(' ')
        species_list = speciesNstrain_split[:2]
        species = ' '.join(species_list)
        strain_list = speciesNstrain_split[2:]
        strain = ' '.join(strain_list)

    return strain.strip()


def save_information_by_tax_level_from_file(taxonomy_level: str, path: str):

    d = {'Num_genes': 'Mean_num_genes', 'Num_variant_genes': 'Mean_num_variant_genes',
         'Num_genes_+strand': 'Mean_num_genes_+strand', 'Num_genes_-strand': 'Mean_num_genes_-strand'}

    categories_ordered = list(categories_order.keys())
    position = categories_order[taxonomy_level] + 1
    taxonomy = categories_ordered[:position]
    cap_taxonomy = [level.capitalize() for level in taxonomy]

    strain_path = f'{path}/strain_all.csv'

    with open(strain_path, 'r', encoding="utf8") as species_file:
        # with open(strain_path, 'r') as species_file:
        #species_info = pd.read_excel(species_file)
        species_info = pd.read_csv(species_file, sep=',', decimal=".")

        # if 'Species' not in species_info.columns and 'Species' not in species_info.columns:
        species_info['Species'] = species_info.apply(
            lambda row: separate_species(row['Species&strain']), axis=1)
        species_info['Strain'] = species_info.apply(
            lambda row: separate_strain(row['Species&strain']), axis=1)

        ordered_columns = ["Taxonomy_Id", "Main_genome_description", "Genomes_mean_size", "Genes_mean_size", "Num_genes", "Num_unique_genes",
                           "Num_genes_+strand", "Num_genes_-strand", "Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain", "Genomes"]
        species_info = species_info[ordered_columns]

        changed_columns = ["Taxonomy_Id", "Main_genome_description", "Genomes_mean_size", "Genes_mean_size", "Num_genes", "Num_variant_genes",
                           "Num_genes_+strand", "Num_genes_-strand", "Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain", "Genomes"]
        species_info.columns = changed_columns

        #species_info = pd.read_csv(species_file, sep=';', decimal=",")
        species_info.dropna(inplace=True)
        species_info_level = species_info.groupby(cap_taxonomy).agg(
            {'Genomes_mean_size': 'mean', 'Genes_mean_size': 'mean', 'Num_genes': 'mean', 'Num_variant_genes': 'mean', 'Num_genes_+strand': 'mean', 'Num_genes_-strand': 'mean', 'Genomes': concat_genomes}).rename(columns=d)

        species_info_level['Num_genomes'] = np.where(species_info_level['Genomes'].str.contains('__'), 0, species_info_level['Genomes'].str.count(
            '-')+1)
        species_info_level.reset_index(inplace=True)
        species_info_level.to_excel(
            f"{path}/{taxonomy_level}_all.xlsx", sheet_name=f'{taxonomy_level}', index=False)

    species_info.to_excel(
        f"{path}/strain_all_separated.xlsx", sheet_name=f'strain', index=False)


args = execute_parser()
path = args.directory

taxonomy_levels = ["superkingdom", "phylum",
                   "class", "order", "family", "genus", "species"]

for taxonomy_level in taxonomy_levels:
    save_information_by_tax_level_from_file(taxonomy_level, path)
