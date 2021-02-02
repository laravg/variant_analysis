import os
import time
from pathlib import Path
from datetime import datetime
import traceback
from Bio import SeqIO
import logging
from modules import input_parser, logic, downloads, outputs, gene_statistics, file_manager
from classes.Genome import Genome, GenomeBatch

taxid_genbankid = {}
genes16s = {}
genes16s_strand = {}
genes16s_begin_end = {}
genes_variants = {}
taxonomy = {}
output_files_ids = {}
total_number_genes16s = {}
unique_number_genes16s = {}
variants_gb_id = {}
variants_strand = {}
id_list = []  # lista de GBID que sube la persona
rejected_ids = []  # GBID que no son bacterias o archaea
failed_ids = []
statistics1_data = {}
statistics2_organism_descr = {}  # Con taxID
statistics2_genome_size = {}  # Con GBID


categories = ["superkingdom", "phylum", "class", "order",
              "family", "genus", "species&strain", "variant"]

ncbi_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id={}&retmode=json&api_key={}&email={}&tool={}"


def main():

    logging.basicConfig(filename='example.log', level=logging.INFO,
                        format='%(levelname)s - %(asctime)s: %(message)s', datefmt='%d/%m/%Y %H:%M:%S %p')
    args = input_parser.execute_parser()

    input_path = args.input
    output_path = args.output
    output_formats = args.format
    training_db = args.training_database

    file_manager.create_dir(output_path)

    result_files = ['genomes', 'genes',
                    'statistics1', 'statistics2'] + output_formats

    for file in result_files:
        dir_path = f'{output_path}/{file}'
        file_manager.create_dir(dir_path)

    if all(x in output_formats for x in ['mothur', 'qiime']):
        dir_path = f'{output_path}/fasta'
        file_manager.create_dir(dir_path)

    logging.info('Starting analysis')

    logging.info('Reading input')

    genomes = []

    if os.path.isfile(input_path):
        genome_ids = file_manager.read_input_file(input_path)
        genomes = downloads.get_genomes(genome_ids)

    elif os.path.isdir(input_path):
        genomes = file_manager.read_input_dir(input_path)

    logging.info('Obtaining genome info from NCBI')
    # Maybe add this methods to Genome class?
    genomes = downloads.get_info_genomes(genomes)
    genome_batches = logic.group_genomes_by_taxonomy(genomes)

    # if rejected_ids:
    #     print(
    #         f'Invalid ids found (not Bacteria or Archaea): {", ".join(rejected_ids)}')
    #     exit(-1)

    logging.info('Analyzing genomes grouped by taxonomy')
    for batch in genome_batches:
        genomes = batch.genomes
        taxonomy_variants = {}
        genes16s_begin_end = {}

        # se obtiene la taxonomía primero para asegurar que no se intenta

        # valid_ids = [
        #     i for i in id_list if i not in rejected_ids]
        # obtener el genoma o los genes 16s de algo que no sea bacteria
        # si hay al menos un identificador válido

        logging.info(f'{batch.genomes_id_str} - Obtaining 16SRNA genes')
        batch = logic.find_genes16s_in_batch(batch, training_db, output_path)

        logging.info(f'{batch.genomes_id_str} - Obtaining variants')
        batch = logic.genes16s_comparison(batch)

        logging.info(
            f'{batch.genomes_id_str} - Generating output files')
        # se generan los identificadores para los archivos de salida
        outputs.generate_id(taxonomy_variants,
                            genes_variants, output_files_ids)
        outputs.save_output(taxonomy_variants, output_files_ids, genes_variants,
                            output_formats, taxid_genbankid, output_path)  # se generan los ficheros salida

        logging.info(f'{batch.genomes_id_str} - Generating statistics')
        gene_statistics.genes16s_statistics(genomes, taxid_genbankid, genes16s, genes_variants, output_files_ids, taxonomy_variants,
                                            variants_gb_id, categories, unique_number_genes16s, total_number_genes16s,
                                            genes16s_strand, genes16s_begin_end, output_path,
                                            statistics1_data, statistics2_organism_descr, statistics2_genome_size)  # se generan los ficheros de estadísticas

        logging.info(f'{batch.genomes_id_str}  analysis completed!')

        # except Exception as e:
        #     print('Fallo ' + str(e))
        #     traceback.print_exception(type(e), e, e.__traceback__)
        #     with open('{}/failed.csv'.format(output_path), 'a+') as failed_ids_path:
        #         failed_ids_path.write('{}\n'.format(genomes_id_str))
        #     with open('{}/log.txt'.format(output_path), 'a+') as log:
        #         log.write(
        #             '{} analysis failed during step \'{}\', continuing\n'.format(genomes_id_str, step))
        #     print("\n\n{} analysis failed, continuing".format(genomes_id_str))

    logging.info('Analysis finished')


main()
