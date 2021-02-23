import os
import sys
import logging
from modules import input_parser, logic, downloads, output_manager, file_manager


def main():

    file_handler = logging.FileHandler(filename='tmp.log')
    file_handler.setFormatter(logging.Formatter(
        '%(levelname)s - %(asctime)s - %(message)s', datefmt='%d/%m/%Y %H:%M:%S'))
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(levelname)s: %(message)s', datefmt='%H:%M:%S'))
    handlers = [file_handler, stdout_handler]
    logging.basicConfig(level=logging.INFO, handlers=handlers)

    args = input_parser.execute_parser()

    input_path = args.input
    output_path = args.output
    output_formats = args.format
    training_db = args.training_database
    taxonomy_level = args.taxonomy_level

    logging.info('Starting analysis')

    logging.info('Reading genomes and downloading from NCBI if needed')

    genomes = []

    try:
        if os.path.isfile(input_path):
            genome_ids, msg = file_manager.read_input_file(input_path)
            if msg:
                logging.warning(msg)
            genomes = downloads.get_genomes(genome_ids)

        elif os.path.isdir(input_path):
            genomes, msg = file_manager.read_input_dir(input_path)
            if msg:
                logging.warning(msg)

    except Exception as e:
        # if no valid genomes
        logging.error(str(e))
        sys.exit()

    result_dirs = ['genes', 'variants_info', 'genes_info'] + output_formats
    file_manager.create_dir(output_path)

    for dir in result_dirs:
        dir_path = f'{output_path}/{dir}'
        file_manager.create_dir(dir_path)

    if all(x in output_formats for x in ['mothur', 'qiime']):
        dir_path = f'{output_path}/fasta'
        file_manager.create_dir(dir_path)

    logging.info('Obtaining genome info from NCBI')

    genomes, msg = downloads.get_info_genomes(genomes)
    if msg:
        logging.warning(msg)
    genome_batches = logic.group_genomes_by_taxonomy(genomes, taxonomy_level)

    logging.info('Analyzing genomes grouped by taxonomy')
    for batch in genome_batches:
        genomes = batch.genomes

        logging.info(f'{batch.genomes_id_str} - Obtaining 16SRNA genes')
        batch = logic.find_genes16s_in_batch(batch, training_db, output_path)

        logging.info(f'{batch.genomes_id_str} - Obtaining variants')
        batch = logic.find_variants(batch)

        logging.info(
            f'{batch.genomes_id_str} - Generating output files')
        output_manager.save_variants(batch, output_formats, output_path)

        # TODO: generate the info for different taxonomic levels
        logging.info(
            f'{batch.genomes_id_str} - Generating file with variants relevant information')
        output_manager.save_variants_information(batch, output_path)
        logging.info(
            f'{batch.genomes_id_str} - Generating file with genes relevant information')
        output_manager.save_genes_information(batch, output_path)

        logging.info(f'{batch.genomes_id_str} - analysis completed')

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
