from classes.Genome import GenomeBatch
from modules.search_16S import initial_motif, final_motif


categories = ["Superkingdom", "Phylum", "Class", "Order",
              "Family", "Genus", "Species&strain", "Variant"]


def save_variants(genome_batch: GenomeBatch, output_formats: list[str], output_dir_path: str):

    output_format_info = {'mothur': 'taxonomy',
                          'dada2': 'fa', 'qiime': 'taxonomy'}
    qiime_categories = ["k__", "p__",
                        "c__", "o__", "f__", "g__", "s__"]

    for format in output_formats:
        path = f'{output_dir_path}/{format}/{genome_batch.genomes_id_str}.{output_format_info[format]}'

        if format == 'dada2':
            with open(path, 'w+') as taxonomy_file:
                lines = ''
                for variant_seq in genome_batch.variants:
                    variant = genome_batch.variants[variant_seq][0]
                    lines += f'>{variant.variant_id}\t{";".join(variant.updated_taxonomy)}\n{variant.sequence}\n'
                taxonomy_file.write(lines)

        else:

            with open(path, 'w+') as taxonomy_file:
                lines = ''
                for variant_seq in genome_batch.variants:
                    variant = genome_batch.variants[variant_seq][0]
                    if format == 'qiime':
                        taxonomy = variant.updated_taxonomy.copy()[:-1]
                        taxonomy[-1] = f'{taxonomy[-1]} {variant.variant_tax_level}'
                        taxonomy = [f'{category}{taxonomy[num]}' for num, category in enumerate(
                            qiime_categories)]
                    else:
                        taxonomy = variant.updated_taxonomy
                    lines += f'>{variant.variant_id}\t{";".join(taxonomy)}\n'
                taxonomy_file.write(lines)

            fasta_dir_path = f'{output_dir_path}/fasta/{genome_batch.genomes_id_str}.fa'
            with open(fasta_dir_path, 'w+') as fasta_file:
                lines = ''
                for variant_seq in genome_batch.variants:
                    variant = genome_batch.variants[variant_seq][0]
                    lines += f'>{variant.variant_id}|{variant.updated_taxonomy[-2]}|{variant.updated_taxonomy[-1]}\n{variant.sequence}'
                fasta_file.write(lines)


def save_variants_information(genome_batch: GenomeBatch, output_dir_path: str):

    header = ["Taxonomy_Id", "New_Id", "Genomes", "Main_Genome_Description",
              "Num_genes", "Num_unique_genes"] + categories + ["Bacteria/Archaea"]

    variants_info = f'{output_dir_path}/variants_info/{genome_batch.genomes_id_str}.csv'

    with open(variants_info, 'w+') as statistics_file:
        statistics_file.write(f'{";".join(header)}\n')
        for variant_seq in genome_batch.variants:  # busca por taxID
            line = []
            variant = genome_batch.variants[variant_seq][0]
            line.append(genome_batch.taxonomy_identifier)
            line.append(variant.variant_id)
            line.append(genome_batch.genomes_id_str)
            line.append(genome_batch.get_main_genome().database_name)
            line.append(str(genome_batch.get_num_analyzed_genes()))
            line.append(str(genome_batch.get_num_variants()))
            line.extend(variant.updated_taxonomy)
            # line = line + tax_variant[:-1]
            # if len(tax_variant[:-1]) < 7:
            #     dash_to_insert = 7 - len(tax_variant[:-1])
            #     line = line + ['_']*dash_to_insert
            # line.append(tax_variant[-1])
            line.append(str(
                genome_batch.taxonomy[0] == 'Bacteria' or genome_batch.taxonomy[0] == 'Archaea'))
            statistics_file.write(f'{";".join(line)}\n')


def save_genes_information(genome_batch: GenomeBatch, output_dir_path: str):

    header = ["Taxonomy ID", "New ID", "GB ID", "GB ID Description", "Genome size", "Initial motif", "Final motif", "Strand +/-", "Variant No.",
              "Initial position", "Final position", "Variant size"] + categories + ["Bacteria/Archaea"]

    statistics2_path = f'{output_dir_path}/genes_info/{genome_batch.genomes_id_str}.csv'

    with open(statistics2_path, 'w+') as statistics_file:
        statistics_file.write(f'{";".join(header)}\n')
        if genome_batch.analyzed_genes:
            for gene in genome_batch.analyzed_genes:
                line = [genome_batch.taxonomy_identifier]
                line.append(gene.variant_id)
                line.append(genome_batch.genomes_id_str)
                line.append(genome_batch.get_main_genome().database_name)
                line.append(
                    str(genome_batch.get_main_genome().database_length))
                line.append(initial_motif)
                line.append(final_motif)
                line.append(gene.strand)
                line.append(gene.variant_tax_level)
                line.append(str(gene.init_position))
                line.append(str(gene.end_position))
                line.append(str(gene.get_seq_length()))
                # line = line + genome.taxonomy
                # if len(genome.taxonomy) < 7:
                #     dash_to_insert = 7 - len(genome.taxonomy)
                #     line = line + ['_']*dash_to_insert
                line.extend(gene.updated_taxonomy)
                line.append(str(
                    genome_batch.taxonomy[0] == 'Bacteria' or genome_batch.taxonomy[0] == 'Archaea'))
                statistics_file.write(f'{";".join(line)}\n')

        else:
            line = [genome_batch.taxonomy_identifier,
                    '_', genome_batch.genomes_id_str]
            line.append(genome_batch.get_main_genome().database_name)
            line.append(str(genome_batch.get_main_genome().database_length))
            line.append(initial_motif)
            line.append(final_motif)
            line.append('_')
            line.append('_')
            line.append('_')
            line.append('_')
            line.append('_')
            line.extend(genome_batch.taxonomy)
            line.append('_')
            line.append(str(
                genome_batch.taxonomy[0] == 'Bacteria' or genome_batch.taxonomy[0] == 'Archaea'))
            statistics_file.write(f'{";".join(line)}\n')
