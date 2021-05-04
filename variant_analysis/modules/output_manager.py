from classes.Genome import GenomeBatch, Genome
from modules.search_16S import initial_motif, final_motif
from modules.logic import categories_order
import numpy as np
import pandas as pd
import os

categories = ["Superkingdom", "Phylum", "Class", "Order",
              "Family", "Genus", "Species&strain", "Variant"]


def save_repeated_genomes(repeated_genomes: list[list[str]], path: str):
    lines = ''
    for genomes in repeated_genomes:
        lines = f'{lines}{", ".join(genomes)}\n'

    with open(path, 'w+') as repeated_genomes_file:
        repeated_genomes_file.write(lines)


def save_repeated_variants(repeated_variants: list[list[str]], path: str):
    lines = ''
    for seq, batches in repeated_variants:
        batches_ids = [batch.genomes_id_str for batch in batches]
        batches_tax = [batch.taxonomy_str for batch in batches]
        lines = f'{lines}{", ".join(batches_tax)} {", ".join(batches_ids)}\n{seq}\n'

    with open(path, 'w+') as repeated_variants_file:
        repeated_variants_file.write(lines)


def save_genomes_without_genes(genomes_without_genes: list[Genome], no_genes_path: str):
    lines = ''
    for genome in genomes_without_genes:
        lines = f'{lines}{genome.identifier}\n'

    with open(no_genes_path, 'w+') as repeated_genomes_file:
        repeated_genomes_file.write(lines)


def save_variants(genome_batch: GenomeBatch, output_formats: list[str], output_dir_path: str):

    output_format_info = {'mothur': 'taxonomy',
                          'dada2': 'fa', 'qiime': 'taxonomy'}
    qiime_categories = {"superkingdom": "k__", "phylum": "p__", "class": "c__",
                        "order": "o__", "family": "f__", "genus": "g__", "species": "s__"}

    for format in output_formats:
        path = f'{output_dir_path}/{format}/{genome_batch.genomes_id_str}.{output_format_info[format]}'

        if format == 'dada2':
            with open(path, 'w+') as taxonomy_file:
                lines = ''
                for variant_seq in genome_batch.variants:
                    variant = genome_batch.variants[variant_seq][0]
                    lines += f'>{variant.variant_id}\t{variant.get_updated_taxonomy_str()}\n{variant.sequence}\n'
                taxonomy_file.write(lines)

        else:

            with open(path, 'w+') as taxonomy_file:
                lines = ''
                for variant_seq in genome_batch.variants:
                    variant = genome_batch.variants[variant_seq][0]
                    if format == 'qiime':
                        taxonomy = variant.updated_taxonomy.copy()
                        taxonomy.pop('variant')
                        taxonomy['species'] = f'{taxonomy["species"]} {variant.variant_tax_level}'
                        taxonomy = [
                            f'{qiime_categories[level]}{taxonomy[level]}' for level in taxonomy]
                    else:
                        taxonomy = variant.get_updated_taxonomy_str()
                    lines += f'>{variant.variant_id}\t{taxonomy}\n'
                taxonomy_file.write(lines)

            fasta_dir_path = f'{output_dir_path}/fasta/{genome_batch.genomes_id_str}.fa'
            with open(fasta_dir_path, 'w+') as fasta_file:
                lines = ''
                for variant_seq in genome_batch.variants:
                    variant = genome_batch.variants[variant_seq][0]
                    lines += f'>{variant.variant_id}|{variant.updated_taxonomy["species"]}|{variant.updated_taxonomy["variant"]}\n{variant.sequence}\n'
                fasta_file.write(lines)


# Information about each gene after grouping by species (1 line per gene)
def save_genes_information(genome_batch: GenomeBatch, output_dir_path: str):

    header = ["Taxonomy ID", "New ID", "GB ID", "GB ID Description",
              "Genome size", "Strand +/-", "Variant No.", "Variant size"] + categories

    statistics2_path = f'{output_dir_path}/genes_info/{genome_batch.genomes_id_str}.csv'

    with open(statistics2_path, 'w+') as statistics_file:
        statistics_file.write(f'{";".join(header)}\n')
        if genome_batch.analyzed_genes:
            main_genome = genome_batch.get_main_genome()
            for gene in genome_batch.analyzed_genes:
                line = [genome_batch.taxonomy_identifier]
                line.append(gene.variant_id)
                line.append(genome_batch.genomes_id_str)
                line.append(main_genome.database_name)
                line.append(main_genome.database_length)
                line.append(gene.strand)
                line.append(gene.variant_tax_level)
                line.append(gene.get_seq_length())
                # line = line + genome.taxonomy
                # if len(genome.taxonomy) < 7:
                #     dash_to_insert = 7 - len(genome.taxonomy)
                #     line = line + ['_']*dash_to_insert
                line.extend(gene.get_updated_taxonomy_list())
                line = list(map(str, line))
                statistics_file.write(f'{";".join(line)}\n')

        else:
            main_genome = genome_batch.get_main_genome()
            line = [genome_batch.taxonomy_identifier, '_', genome_batch.genomes_id_str, main_genome.database_name,
                    main_genome.database_length, '_', '_', '_', genome_batch.get_taxonomy_list(), '_']
            line = list(map(str, line))
            statistics_file.write(f'{";".join(line)}\n')


# Information about each variant after grouping by species (1 line per variant)
def save_variants_information(genome_batch: GenomeBatch, output_dir_path: str):

    header = ["Taxonomy_Id", "New_Id", "Genomes",
              "Main_genome_description", "Num_genes", "Num_copies"] + categories

    variants_info = f'{output_dir_path}/variants_info/{genome_batch.genomes_id_str}.csv'

    with open(variants_info, 'w+') as statistics_file:
        statistics_file.write(f'{";".join(header)}\n')
        num_copies = genome_batch.get_num_copies_per_variant()
        for variant_seq in genome_batch.variants:  # busca por taxID
            line = []
            variant = genome_batch.variants[variant_seq][0]
            line.append(genome_batch.taxonomy_identifier)
            line.append(variant.variant_id)
            line.append(genome_batch.genomes_id_str)
            line.append(genome_batch.get_main_genome().database_name)
            line.append(genome_batch.get_num_analyzed_genes())
            line.append(num_copies[variant.updated_taxonomy['variant']])
            line.extend(variant.get_updated_taxonomy_list())
            # line = line + tax_variant[:-1]
            # if len(tax_variant[:-1]) < 7:
            #     dash_to_insert = 7 - len(tax_variant[:-1])
            #     line = line + ['_']*dash_to_insert
            # line.append(tax_variant[-1])
            line = list(map(str, line))
            statistics_file.write(f'{";".join(line)}\n')


# Information about genes after grouping by species (1 line per file), such as the mean number of genes in all genomes analyzed for that species.
def save_species_information(genome_batch: GenomeBatch, output_dir_path: str):

    header = ["Taxonomy_Id", "Main_genome_description", "Genomes_mean_size",
              "Genes_mean_size", "Num_genes", "Num_variant_genes", "Num_genes_+strand", "Num_genes_-strand"] + categories[:-1] + ["Genomes"]

    variants_info = f'{output_dir_path}/7_species_info/{genome_batch.genomes_id_str}.csv'
    all_lines = []
    with open(variants_info, 'w+') as statistics_file:
        if genome_batch.analyzed_genes:
            statistics_file.write(f'{";".join(header)}\n')
            line = []
            line.append(genome_batch.taxonomy_identifier)
            line.append(genome_batch.get_main_genome().database_name)
            line.append(genome_batch.get_genomes_with_genes_mean_size())
            line.append(genome_batch.get_genes_mean_size())
            line.append(genome_batch.get_num_analyzed_genes())
            line.append(genome_batch.get_num_variants())
            line.append(genome_batch.get_num_genes_per_strand('+'))
            line.append(genome_batch.get_num_genes_per_strand('-'))
            line.extend(genome_batch.get_taxonomy_list())
            line.append(genome_batch.get_genomes_with_genes_id_str())
            # line = line + tax_variant[:-1]
            # if len(tax_variant[:-1]) < 7:
            #     dash_to_insert = 7 - len(tax_variant[:-1])
            #     line = line + ['_']*dash_to_insert
            # line.append(tax_variant[-1])
            line = list(map(str, line))
            all_lines.append(line)
            statistics_file.write(f'{";".join(line)}\n')
        else:
            main_genome = genome_batch.get_main_genome()
            line = [genome_batch.taxonomy_identifier, main_genome.database_name,
                    '_', '_', '_', '_', '_', '_'] + genome_batch.get_taxonomy_list() + ['_']
            all_lines.append(line)
            statistics_file.write(f'{";".join(line)}\n')

    all_variants_info_path = f'{output_dir_path}/7_species_info/species_all.csv'
    with open(all_variants_info_path, 'a+') as statistics_file:
        if os.stat(all_variants_info_path).st_size == 0:
            statistics_file.write(f'{";".join(header)}\n')
        for line in all_lines:
            statistics_file.write(f'{";".join(line)}\n')


# Information about genes after grouping by selected tax_level (1 line per file, 1 file per genome batches grouped with selected tax
# level), such as the mean number of genes in all genomes analyzed for that species. If one batch contains genomes without any genes, that genome will be excluded.
def save_information_by_tax_level(grouped_batches: dict[str, list[GenomeBatch]], taxonomy_level: str, output_dir_path: str):

    position = categories_order[taxonomy_level] + 1
    header = ["Tax_level", "Num_genomes", "Genomes_mean_size", "Genes_mean_size", "Mean_num_genes",
              "Mean_num_variant_genes", "Mean_num_genes_+strand", "Mean_num_genes_-strand"] + categories[:position] + ["Genomes"]
    all_lines = []
    for taxonomy, batch_group in grouped_batches.items():
        variants_info_path = f'{output_dir_path}/{categories_order[taxonomy_level]+1}_{taxonomy_level}_info/{taxonomy}.csv'
        batches_with_genes = [
            batch for batch in batch_group if batch.analyzed_genes]
        if batches_with_genes:
            with open(variants_info_path, 'w+') as statistics_file:
                statistics_file.write(f'{";".join(header)}\n')
                line = []
                line.append(taxonomy_level)
                num_genomes = np.sum(
                    [1 for batch in batch_group for genome in batch.genomes if genome.amplified_genes])
                line.append(np.around(num_genomes, decimals=3))
                genomes_mean_size = np.mean([batch.get_genomes_with_genes_mean_size()
                                             for batch in batch_group if batch.analyzed_genes])
                line.append(np.around(genomes_mean_size, decimals=3))
                genes_mean_size = np.mean([batch.get_genes_mean_size()
                                           for batch in batch_group if batch.analyzed_genes])
                line.append(np.around(genes_mean_size, decimals=3))
                num_analyzed_genes = np.mean(
                    [batch.get_num_analyzed_genes() for batch in batch_group if batch.analyzed_genes])
                line.append(np.around(num_analyzed_genes, decimals=3))
                num_variants = np.mean([batch.get_num_variants()
                                        for batch in batch_group if batch.analyzed_genes])
                line.append(np.around(num_variants, decimals=3))
                num_genes_per_strand_pos = np.mean(
                    [batch.get_num_genes_per_strand('+') for batch in batch_group if batch.analyzed_genes])
                line.append(np.around(num_genes_per_strand_pos, decimals=3))
                num_genes_per_strand_neg = np.mean(
                    [batch.get_num_genes_per_strand('-') for batch in batch_group if batch.analyzed_genes])
                line.append(np.around(num_genes_per_strand_neg, decimals=3))
                line.append(taxonomy)
                genome_ids_str = '-'.join(
                    [batch.get_genomes_with_genes_id_str() for batch in batch_group if batch.analyzed_genes])
                line.append(genome_ids_str)

                line = list(map(str, line))
                all_lines.append(line)
                statistics_file.write(f'{";".join(line)}\n')

    all_variants_info_path = f'{output_dir_path}/{categories_order[taxonomy_level]+1}_{taxonomy_level}_info/{taxonomy_level}_all.csv'
    with open(all_variants_info_path, 'w+') as statistics_file:
        statistics_file.write(f'{";".join(header)}\n')
        for line in all_lines:
            statistics_file.write(f'{";".join(line)}\n')


# Information about genes after grouping by selected tax_level (1 line per file, 1 file per genome batches grouped with selected tax
# level), such as the mean number of genes in all genomes analyzed for that species. If one batch contains genomes without any genes, that genome will be excluded.
def save_information_by_tax_level_from_file(taxonomy_level: str, species_path: str, output_dir_path: str):
    d = {'Num_genes': 'Mean_num_genes', 'Num_variant_genes': 'Mean_num_variant_genes',
         'Num_genes_+strand': 'Mean_num_genes_+strand', 'Num_genes_-strand': 'Mean_num_genes_-strand'}

    categories_ordered = list(categories_order.keys())
    position = categories_order[taxonomy_level] + 1
    taxonomy = categories_ordered[:position]
    cap_taxonomy = [level.capitalize() for level in taxonomy]

    with open(species_path, 'r', encoding="utf8") as species_file:
        species_info = pd.read_csv(species_file, sep=';', decimal=",")
        species_info.dropna(inplace=True)
        species_info_level = species_info.groupby(cap_taxonomy).agg(
            {'Genomes_mean_size': 'mean', 'Genes_mean_size': 'mean', 'Num_genes': 'mean', 'Num_variant_genes': 'mean', 'Num_genes_+strand': 'mean', 'Num_genes_-strand': 'mean', 'Genomes': '-'.join}).rename(columns=d)
        species_info_level['Num_genomes'] = species_info_level['Genomes'].str.count(
            '-')+1
        species_info_level.reset_index(inplace=True)
        species_info_level.to_excel(
            f"resultsfileBACTERIA/summary/rerun/{taxonomy_level}_all_rerun.xlsx", sheet_name=f'{taxonomy_level}', index=False)
