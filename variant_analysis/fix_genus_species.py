import pandas as pd
import numpy as np
import argparse


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


args = execute_parser()
path = args.directory

species_df = pd.read_csv(path, sep=',', decimal=".")


species_df['Genus != Species ?'] = np.where(species_df['Genus'].apply(lambda x: x.split(' ')[0]) != species_df['Species&strain'].apply(
    lambda x: x.split(' ')[0]), 'Distinct!', '')

species_df['Species&strainFixed'] = np.where(species_df['Genus != Species ?'].str.contains("Distinct!"), species_df['Genus'] + ' ' +
                                             species_df['Species&strain'].apply(lambda x: ' '.join(x.split(' ')[1:])), species_df['Species&strain'])

new_names = ["Taxonomy_Id", "Main_genome_description", "Genomes_mean_size", "Genes_mean_size", "Num_genes", "Num_unique_genes", "Num_genes_+strand",
             "Num_genes_-strand", "Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species&strain", 'Genus != Species ?', 'Species&strainFixed', "Genomes"]
species_df = species_df[new_names]


new_path = path.replace('.csv', '')
with pd.ExcelWriter(f'{new_path}_fix_genus_species.xlsx', engine='xlsxwriter') as writer:
    species_df.to_excel(writer, sheet_name='species_all',
                        index=None, header=True)

    # Get the xlsxwriter workbook and worksheet objects.
    workbook = writer.book
    worksheet = writer.sheets['species_all']

    # Get the dimensions of the dataframe.
    (max_row, max_col) = species_df.shape

    # Create a list of column headers, to use in add_table().
    column_settings = []
    for header in species_df.columns:
        column_settings.append({'header': header})

    # Add the table.
    worksheet.add_table(0, 0, max_row, max_col - 1,
                        {'columns': column_settings})

    # Make the columns wider for clarity.
    worksheet.set_column(9, 13, 18)
    worksheet.set_column(14, 14, 45)
    worksheet.set_column(16, 16, 45)
