import pandas as pd
import numpy as np
import argparse
import os
from pathlib import Path


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

categories = ["superkingdom_all.xlsx", "phylum_all.xlsx", "class_all.xlsx", "order_all.xlsx",
              "family_all.xlsx", "genus_all.xlsx", "species_all.xlsx", "strain_all.xlsx"]
# "species_all_fix_genus_species.xlsx"


summary_path = Path(path)
new_file_path = summary_path / 'all.xlsx'
writer = pd.ExcelWriter(new_file_path, engine='xlsxwriter')

for file in categories:

    file_path = summary_path / file
    df = pd.read_excel(file_path)
    #df.rename(columns={'Num_unique_genes': 'Num_variant_genes'}, inplace=True)

    # if file == 'species_all_fix_genus_species.xlsx':
    #     sheet_name = file.replace('_all_fix_genus_species.xlsx', '').title()
    #     columns = ["Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species&strain", "Genomes_mean_size", "Genes_mean_size",
    #                "Num_genes", "Num_variant_genes", "Num_genes_+strand", "Num_genes_-strand", "Genomes", "Main_genome_description", "Taxonomy_Id"]
    # else:
    #     sheet_name = file.replace('_all.xlsx', '').title()
    #     columns = list(df.columns)
    #     columns[-1], columns[-2] = columns[-2], columns[-1]

    sheet_name = file.replace('_all.xlsx', '').title()
    columns = list(df.columns)
    columns[-1], columns[-2] = columns[-2], columns[-1]

    new_df = df[columns]
    new_df.to_excel(writer, sheet_name=sheet_name,
                    index=None, header=True, columns=columns)
    workbook = writer.book
    worksheet = writer.sheets[sheet_name]

    # Get the dimensions of the dataframe.
    (max_row, max_col) = new_df.shape

    # Create a list of column headers, to use in add_table().
    column_settings = []
    for header in new_df.columns:
        column_settings.append({'header': header})

    # Add the table.
    worksheet.add_table(0, 0, max_row, max_col - 1,
                        {'columns': column_settings})

    # Make the columns wider for clarity.
    worksheet.set_column(0, max_col - 1, 20)

    format_align = workbook.add_format()
    format_align.set_align('center')
    format_align.set_align('vcenter')

    worksheet.set_column(0, max_col - 2, 20, format_align)

    format_decimal = workbook.add_format()
    format_decimal.set_num_format('0.000')
    format_decimal.set_align('center')
    format_decimal.set_align('vcenter')

    # if file == 'species_all_fix_genus_species.xlsx':
    if file == 'strain_all.xlsx':
        worksheet.set_column(2, 3, 20, format_decimal)
    else:
        worksheet.set_column(max_col - 8, max_col - 3, 20, format_decimal)
        worksheet.set_column(max_col - 1, max_col - 1, 100)


writer.save()
writer.close()
