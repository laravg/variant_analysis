import pathlib
import re
import os
import urllib.request
import modules.config as config
import time
from datetime import datetime
from datetime import datetime
from modules.search_16S import initial_motif, final_motif
from bs4 import BeautifulSoup
import json

ncbi_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id={}&retmode=json&api_key={}&email={}&tool={}"


def genes16s_statistics(genomes, taxid_genbankid, genes16s, genes_variants, output_files_ids, taxonomy_variants,
                        variants_gb_id, categories, unique_number_genes16s, total_number_genes16s,
                        genes16s_strand, genes16s_begin_end, output_dir_path,
                        statistics1_data, statistics2_organism_descr, statistics2_genome_size):  # se generar los ficheros de estadísticas

    header = ["Taxonomy ID", "New ID", "GB ID", "GB ID Description",
              "No. of genes", "No. of unique genes"] + categories + ["Bacteria/Archaea"]
    header_line = ""
    for item in header:
        if header.index(item) == 0:
            header_line = "{}".format(item)
        else:
            header_line = header_line + ";{}".format(item)
    # primera tabla, se muestra información sobre las variantes encontradas
    GBIDs = '_'.join([x.identifier for x in genomes])
    statistics1_path = '{}/statistics1/{}.csv'.format(output_dir_path, GBIDs)

    with open(statistics1_path, 'w+') as statistics_file:
        statistics_file.write(header_line + '\n')
        line = []
        for item in taxonomy_variants:  # busca por taxID
            for tax_variant in taxonomy_variants.get(item):
                line.append(item)
                line.append(output_files_ids.get(item)[tax_variant[-1]])
                line.append(variants_gb_id[item])
                line.append(statistics1_data[item])
                line.append(total_number_genes16s.get(item))
                line.append(unique_number_genes16s.get(item))
                line = line + tax_variant[:-1]
                if len(tax_variant[:-1]) < 7:
                    dash_to_insert = 7 - len(tax_variant[:-1])
                    line = line + ['_']*dash_to_insert
                line.append(tax_variant[-1])
                if tax_variant[0] == 'bacteria' or 'archaea':
                    line.append('Si')
                else:
                    line.append('No')
                entry = ''
                for field in line:
                    entry = '{0}{1};'.format(entry, field)
                entry = entry + '\n'
                statistics_file.write(entry)
                line = []

    # segunda tabla: se muestran información sobre todos identificadores indicados y sus genes si los tiene
    header = ["Taxonomy ID", "New ID", "GB ID", "GB ID Description", "Genome size", "Initial motif", "Final motif", "Strand +/-", "Variant No.",
              "Initial position", "Final position", "Variant size"] + categories + ["Bacteria/Archaea"]
    header_line = ""
    for item in header:
        if header.index(item) == 0:
            header_line = "{}".format(item)
        else:
            header_line = header_line + ";{}".format(item)

    statistics2_path = '{}/statistics2/{}.csv'.format(output_dir_path, GBIDs)

    with open(statistics2_path, 'w+') as statistics_file:
        statistics_file.write(header_line + '\n')

        for genome in genomes:  # obtiene GBID que subio la persona
            GBID = genome.identifier
            taxonomy_id = genome.taxid

            if GBID in genes16s:  # si el GBID tiene genes 16s y es Bacteria o Archaea
                gene_position = 0
                # recorre los genes 16s de ese GBID
                for gene16s in genes16s[GBID]:
                    line = [taxonomy_id]
                    variant = None
                    # obtiene las variantes de ese gen16s
                    for var in genes_variants[taxonomy_id]:
                        if genes_variants[taxonomy_id][var] == gene16s:
                            # se anhade la columna de id del programa
                            line.append(output_files_ids.get(taxonomy_id)[var])
                            variant = var
                    line.append(GBID)  # se anhade la columna de GBID
                    # se anhade la columna de descripcion de GBID
                    line.append(statistics2_organism_descr[GBID])
                    # se anhade la columna de tamanho de genoma
                    line.append(statistics2_genome_size[GBID])
                    # se anhade la columna de motif inicial
                    line.append(initial_motif)
                    # se anhade la columna de motif final
                    line.append(final_motif)
                    # se anhade la columna de hebra
                    line.append(genes16s_strand[GBID][gene_position])
                    line.append(variant)  # se anhade la columna de variante
                    # se anhade la columna de initial position
                    line.append(
                        genes16s_begin_end[GBID][gene_position][0])
                    # se anhade la columna de final position
                    line.append(
                        genes16s_begin_end[GBID][gene_position][1])
                    # se anhade la columna de tamanho de la variante
                    line.append(len(gene16s))
                    # se anhaden las columnas de taxonomia + variante
                    line = line + genome.taxonomy
                    if len(genome.taxonomy) < 7:
                        dash_to_insert = 7 - len(genome.taxonomy)
                        line = line + ['_']*dash_to_insert
                    line.append(variant)
                    # se anhade la columna de Bacteria/Archaea
                    line.append('Yes')
                    entry = ''

                    for field in line:
                        entry = '{0}{1};'.format(entry, field)
                    entry = entry + '\n'
                    statistics_file.write(entry)
                    gene_position = gene_position + 1

            # hacer fila en caso de que no tenga genes16s
            elif GBID not in genes16s:
                # se anhaden las columnas de taxID, id del programa y GBID
                line = [taxonomy_id, '_', GBID]
                # se anhade la columna de descripcion de GBID
                line.append(statistics2_organism_descr[GBID])
                # se anhade la columna de tamanho de genoma
                line.append(statistics2_genome_size[GBID])
                # se anhade la columna de motif inicial
                line.append(initial_motif)
                line.append(final_motif)  # se anhade la columna de motif final
                line.append('_')  # se anhade la columna de hebra
                line.append('_')  # se anhade la columna de variante
                line.append('_')  # se anhade la columna de initial position
                line.append('_')  # se anhade la columna de final position
                # se anhade la columna de tamanho de la variante
                line.append('_')
                # se anhaden las columnas de taxonomia
                line = line + genome.taxonomy
                if len(genome.taxonomy) < 7:
                    dash_to_insert = 7 - len(genome.taxonomy)
                    line = line + ['_']*dash_to_insert
                line.append('_')
                line.append('Yes')  # se anhade la columna de Bacteria/Archaea
                entry = ''

                for field in line:
                    entry = '{0}{1};'.format(entry, field)
                entry = entry + '\n'
                statistics_file.write(entry)

            # hacer fila en caso de que no sea bacteria o archaea
            # elif GBID in rejected_ids:
            #     # se anhaden las columnas de taxID, id del programa y GBID
            #     line = [taxonomy_id, '_', GBID]
            #     # se anhade la columna de descripcion de GBID
            #     line.append(statistics2_organism_descr[GBID])
            #     # se anhade la columna de tamanho de genoma
            #     line.append(statistics2_genome_size[GBID])
            #     # se anhade la columna de motif inicial
            #     line.append(initial_motif)
            #     line.append(final_motif)  # se anhade la columna de motif final
            #     line.append('_')  # se anhade la columna de hebra
            #     line.append('_')  # se anhade la columna de variante
            #     line.append('_')  # se anhade la columna de initial position
            #     line.append('_')  # se anhade la columna de final position
            #     # se anhade la columna de tamanho de la variante
            #     line.append('_')
            #     for i in categories:
            #         line.append('_')  # se anhaden las columnas de taxonomia
            #     line.append('No')  # se anhade la columna de Bacteria/Archaea
            #     entry = ''

            #     for field in line:
            #         entry = '{0}{1};'.format(entry, field)
            #     entry = entry + '\n'
            #     statistics_file.write(entry)
