import os
import sys
import logging
from modules import input_parser, logic, downloads, output_manager, file_manager
from modules.logic import categories_order


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

    input_path: str = args.input
    output_path: str = args.output
    output_formats: list[str] = args.format
    training_db: str = args.training_database
    taxonomy_levels: list[str] = args.taxonomy_level
    substitute_iupac_path: str = args.substitute_iupac
    primers_path: str = args.primers

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

    if substitute_iupac_path:
        file_manager.create_dir(substitute_iupac_path)
        genomes = logic.substitute_IUPAC_genome_sequences(
            genomes, substitute_iupac_path)

    primer_pairs = []
    if primers_path:
        primer_pairs, msg = file_manager.read_csv_file(primers_path)
        if msg:
            logging.warning(msg)

    result_dirs = ['genes', 'variants_info', 'genes_info', '7_species_info'] + \
        [f'{categories_order[taxonomy_level]+1}_{taxonomy_level}_info' for taxonomy_level in taxonomy_levels if taxonomy_level != 'species'] + \
        output_formats

    logging.info('Obtaining genome info from NCBI')

    genomes, msg = downloads.get_info_genomes(genomes)
    if msg:
        logging.warning(msg)

    repeated_genomes = logic.check_repeated_genomes(genomes)

    genome_batches = logic.group_genomes_by_taxonomy(genomes)

    logging.info('Analyzing genomes grouped by taxonomy')

    if primer_pairs:
        root_output_path = output_path
        for primer_pair in primer_pairs:
            logic.reset_genome_batches(genome_batches)

            output_path = f'{root_output_path}_{primer_pair.identifier}'
            file_manager.create_dir(output_path)

            for dir in result_dirs:
                dir_path = f'{output_path}/{dir}'
                file_manager.create_dir(dir_path)

            if all(x in output_formats for x in ['mothur', 'qiime']):
                dir_path = f'{output_path}/fasta'
                file_manager.create_dir(dir_path)

            for batch in genome_batches:
                logging.info(
                    f'{batch.genomes_id_str} - Obtaining 16SRNA genes with primer pair {primer_pair.identifier}')
                batch = logic.trim_genes16s_in_batch(
                    batch, primer_pair, output_path)

                logging.info(f'{batch.genomes_id_str} - Obtaining variants')
                batch = logic.find_variants(batch)

                logging.info(
                    f'{batch.genomes_id_str} - Generating output files')
                output_manager.save_variants(
                    batch, output_formats, output_path)

                logging.info(
                    f'{batch.genomes_id_str} - Generating file with genes relevant information')
                output_manager.save_genes_information(batch, output_path)

                logging.info(
                    f'{batch.genomes_id_str} - Generating file with variants relevant information')
                output_manager.save_variants_information(batch, output_path)

                logging.info(
                    f'{batch.genomes_id_str} - Generating file with species relevant information')
                output_manager.save_species_information(batch, output_path)

                logging.info(f'{batch.genomes_id_str} - analysis completed')

            repeated_variants = logic.check_repeated_variants(genome_batches)
            if repeated_variants:
                path = f'{output_path}/repeated_variants.txt'
                output_manager.save_repeated_variants(repeated_variants, path)
                logging.warning(
                    'Repeated variants by sequence found, file with sequence and identifiers generated')

            genomes_without_genes = logic.check_genomes_without_genes(genomes)
            no_genes_path = f'{output_path}/genomes_without_genes.txt'
            if genomes_without_genes:
                output_manager.save_genomes_without_genes(
                    genomes_without_genes, no_genes_path)

            for taxonomy_level in taxonomy_levels:
                if taxonomy_level != 'species':
                    grouped_batches = logic.group_batches_by_tax_level(
                        genome_batches, taxonomy_level)
                    logging.info(
                        f'Generating file with relevant information about selected taxonomic level: {taxonomy_level}')
                    output_manager.save_information_by_tax_level(
                        grouped_batches, taxonomy_level, output_path)

    else:
        file_manager.create_dir(output_path)
        if repeated_genomes:
            path = f'{output_path}/repeated_genomes.txt'
            output_manager.save_repeated_genomes(repeated_genomes, path)
            logging.error(
                'Repeated genomes by sequence found, file with identifiers generated and analysis continues')

        for dir in result_dirs:
            dir_path = f'{output_path}/{dir}'
            file_manager.create_dir(dir_path)

        if all(x in output_formats for x in ['mothur', 'qiime']):
            dir_path = f'{output_path}/fasta'
            file_manager.create_dir(dir_path)

        for batch in genome_batches:

            logging.info(f'{batch.genomes_id_str} - Obtaining 16SRNA genes')
            batch = logic.find_genes16s_in_batch(
                batch, training_db, output_path)

            logging.info(f'{batch.genomes_id_str} - Obtaining variants')
            batch = logic.find_variants(batch)

            logging.info(
                f'{batch.genomes_id_str} - Generating output files')
            output_manager.save_variants(batch, output_formats, output_path)

            logging.info(
                f'{batch.genomes_id_str} - Generating file with genes relevant information')
            output_manager.save_genes_information(batch, output_path)

            logging.info(
                f'{batch.genomes_id_str} - Generating file with variants relevant information')
            output_manager.save_variants_information(batch, output_path)

            logging.info(
                f'{batch.genomes_id_str} - Generating file with species relevant information')
            output_manager.save_species_information(batch, output_path)

            logging.info(f'{batch.genomes_id_str} - analysis completed')

        repeated_variants = logic.check_repeated_variants(genome_batches)
        if repeated_variants:
            path = f'{output_path}/repeated_variants.txt'
            output_manager.save_repeated_variants(repeated_variants, path)
            logging.warning(
                'Repeated variants by sequence found, file with sequence and identifiers generated')

        genomes_without_genes = logic.check_genomes_without_genes(genomes)
        no_genes_path = f'{output_path}/genomes_without_genes.txt'
        if genomes_without_genes:
            output_manager.save_genomes_without_genes(
                genomes_without_genes, no_genes_path)

        for taxonomy_level in taxonomy_levels:
            if taxonomy_level != 'species':
                grouped_batches = logic.group_batches_by_tax_level(
                    genome_batches, taxonomy_level)
                logging.info(
                    f'Generating file with relevant information about selected taxonomic level: {taxonomy_level}')
                output_manager.save_information_by_tax_level(
                    grouped_batches, taxonomy_level, output_path)

        # Analysis from a existent species file
        # species_path = 'resultsfileBACTERIA/summary/species_all_rerun.csv'
        # for taxonomy_level in taxonomy_levels:
        #     output_manager.save_information_by_tax_level_from_file(taxonomy_level, species_path

        logging.info('Analysis finished')


main()
