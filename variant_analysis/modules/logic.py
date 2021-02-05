from Bio import SeqIO
from modules.search_16S import find
import hashlib
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
        genome_batch.analyzed_genes.extend(genome.amplified_genes)

    return genome_batch


def find_genes16s(genome: Genome, training_db: str, output_dir_path: str) -> Genome:

    genes_dir_path = f'{output_dir_path}/genes/{genome.identifier}.fa'

    kmer_path = f'utilities/kmer{training_db}'
    with open(genes_dir_path, 'w+') as genes:
        find(genome.get_fasta_contents(), kmer_path, genes,
             4, 1000, 2500)  # default values

    with open(genes_dir_path, 'r') as genes:
        fasta_sequences = SeqIO.parse(genes, 'fasta')
        for fasta in fasta_sequences:
            header, sequence = fasta.description, str(
                fasta.seq)
            gene = Gene(header, sequence)
            genome.amplified_genes.append(gene)

    return genome


def find_variants(genome_batch: GenomeBatch) -> GenomeBatch:

    gene_seqs = [gene.sequence for gene in genome_batch.analyzed_genes]
    variant_gene_seqs = set(gene_seqs)
    genome_batch.variants_tax_level = {
        seq: f'v{num}' for num, seq in enumerate(variant_gene_seqs, 1)}

    for gene in genome_batch.analyzed_genes:
        gene.variant_tax_level = genome_batch.variants_tax_level[gene.sequence]
        taxonomy_copy = genome_batch.taxonomy.copy()
        taxonomy_copy.append(gene.variant_tax_level)
        gene.updated_taxonomy = taxonomy_copy
        gene.variant_id = str(int(hashlib.sha256(
            gene.sequence.encode('utf-8')).hexdigest(), 16) % 10**6).zfill(6)

    variants = {}
    for gene in genome_batch.analyzed_genes:
        if gene.sequence not in variants:
            variants[gene.sequence] = []
        variants[gene.sequence].append(gene)

    genome_batch.variants = variants

    return genome_batch
