import json
from Bio import SeqIO, Entrez
from modules.search_16S import find
import os
from modules import config, file_manager
from classes.Genome import Genome, GenomeBatch
from classes.Gene import Gene


def group_genomes_by_taxonomy(genomes: list[Genome]) -> list[GenomeBatch]:

    genomes_grouped = {}
    taxonomy = {}
    for genome in genomes:
        if genome.taxonomy_identifier not in genomes_grouped:
            genomes_grouped[genome.taxonomy_identifier] = []
        genomes_grouped[genome.taxonomy_identifier].append(genome)
        taxonomy[genome.taxonomy_identifier] = genome.taxonomy

    genome_batches = []
    for tax_id in genomes_grouped:
        batch = GenomeBatch(
            tax_id, taxonomy[tax_id], genomes_grouped[tax_id])
        genome_batches.append(batch)

    return genome_batches


def find_genes16s_in_batch(genome_batch: GenomeBatch, training_db: str, output_dir_path: str) -> GenomeBatch:

    for genome in genome_batch.genomes:
        genome = find_genes16s(genome, training_db, output_dir_path)
        genome_batch.genes_to_analyse.extend(genome.amplified_genes)

    return genome_batch


def find_genes16s(genome: Genome, training_db: str, output_dir_path: str) -> Genome:

    genes_dir_path = f'{output_dir_path}/genes/{genome.identifier}.fa'
    kmer_path = f'utilities/kmer{training_db}'

    with open(genes_dir_path, 'w+') as genes:
        find(genome.get_fasta_contents(), kmer_path, genes,
             4, 1000, 2500)  # default values
        fasta_sequences = SeqIO.parse(genes, 'fasta')
        for fasta in fasta_sequences:
            header, sequence = fasta.description, str(
                fasta.seq)
            gene = Gene(header, sequence)
            # if gbid not in valid_ids:  # remove dot
            #     gbid_clean = gbid.split(".")
            #     gbid = gbid_clean[0]
            genome.amplified_genes.append(gene)

    return genome


def genes16s_comparison(genome_batch: GenomeBatch):  # se comparan los genes 16S

    gene_seqs = [gene.sequence for gene in genome_batch.genes_to_analyse]
    variant_gene_seqs = set(gene_seqs)
    genome_batch.variant_gene_seqs_ids = {
        seq: f'v{num}' for num, seq in enumerate(variant_gene_seqs, 1)}

    for gene in genome_batch.genes_to_analyse:
        gene.variant_id = genome_batch.variant_gene_seqs_ids[gene.sequence]
        taxonomy_copy = genome_batch.taxonomy.copy()
        gene.updated_taxonomy = taxonomy_copy.append(gene.variant_id)

    return genome_batch


# def get_gen16s_multiple_files(genome, genes16s, training_db, genes16s_strand, genes16s_begin_end, output_dir_path, file_path):
#     genes_dir_path = '{}/genes/{}'.format(output_dir_path, genome.file_name)
#     genomes_dir_path = '{}/{}'.format(file_path, genome.file_name)
#     kmer_path = "{}/utilities/kmer{}".format(
#         pathlib.Path(__file__).parent.absolute(), training_db)

#     with open(genes_dir_path, 'w+') as genes, open(genomes_dir_path, 'r') as genomes:
#         find(genomes, kmer_path, genes, 4, 1000, 2500)  # default values
#     with open(genes_dir_path, 'r') as genes:
#         lines = genes.read().splitlines()
#         iterator = iter(lines)
#         for line in iterator:
#             split_line = line.split(";")
#             gbid = split_line[0].split(">")
#             gbid = gbid[1]
#             if gbid not in genes16s:
#                 genes16s[gbid] = []
#             if gbid not in genes16s_strand:
#                 genes16s_strand[gbid] = []
#             genes16s[gbid].append(next(iterator))
#             genes16s_strand[gbid].append(split_line[1].split(']')[-1])
#             begin_end = [split_line[1].split(',')[0].strip(
#                 '['), split_line[1].split(',')[1].split(']')[0]]
#             if gbid not in genes16s_begin_end:
#                 genes16s_begin_end[gbid] = []
#             genes16s_begin_end[gbid].append(begin_end)
