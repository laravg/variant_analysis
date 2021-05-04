from classes.PrimerPair import PrimerPair
from Bio import SeqIO
from modules.search_16S import find
import hashlib
from classes.Genome import Genome, GenomeBatch
from classes.Gene import Gene
from classes.PrimerPair import PrimerPair, TargetSequence
from Bio.Data import IUPACData
from Bio.Seq import Seq
import regex as re
import random


categories_order = {"superkingdom": 0, "phylum": 1, "class": 2, "order": 3,
                    "family": 4, "genus": 5, "species": 6}

iupac_code = {'R': ['A', 'G'], 'Y': ['T', 'C'], 'M': ['A', 'C'], 'K': ['T', 'G'], 'S': ['C', 'G'], 'W': ['A', 'T'], 'H': [
    'A', 'C', 'T'], 'B': ['G', 'C', 'T'], 'V': ['A', 'C', 'G'], 'D': ['A', 'G', 'T'], 'N': ['G', 'A', 'C', 'T']}

max_size = 2300
min_size = 100

MISMATCHES = 0


def substituteIUPAC(matchobj):
    return matchobj.group(1) + ''.join([random.choice(iupac_code[char]) for char in matchobj.group(2)]) + matchobj.group(3)


def substitute_IUPAC_genome_sequences(genomes: list[Genome], output_dir_path: str):
    pattern = '([ACGT]|^)([NBVMYSKRWHD]{1,2})([ACGT]|$)'
    for genome in genomes:
        modified_sequence = re.sub(pattern, substituteIUPAC, genome.sequence)
        genome.sequence = modified_sequence

        with open(f'{output_dir_path}/{genome.identifier}.fa', 'w+') as genome_file:
            genome_file.write(f'{genome.get_fasta_contents()}')

    return genomes


def check_repeated_genomes(genomes: list[Genome]):
    sequences = {}
    for genome in genomes:
        hashed_sequence = hashlib.md5(genome.sequence.encode()).hexdigest()
        if hashed_sequence not in sequences:
            sequences[hashed_sequence] = []
        sequences[hashed_sequence].append(genome.identifier)

    repeated_sequences = [ids for ids in list(
        sequences.values()) if len(ids) > 1]

    return repeated_sequences


def check_repeated_variants(genome_batches: list[GenomeBatch]):
    variants = {}
    for batch in genome_batches:
        for variant in list(batch.variants.keys()):
            if variant not in variants:
                variants[variant] = []
            variants[variant].append(batch)

    repeated_variants = [(seq, variants[seq])
                         for seq in variants if len(variants[seq]) > 1]

    return repeated_variants


def check_genomes_without_genes(genomes: list[Genome]):
    genomes_without_genes = [
        genome for genome in genomes if not genome.amplified_genes]
    return genomes_without_genes


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
            gene = Gene(sequence)
            gene.init_from_header(header)
            genome.amplified_genes.append(gene)

    return genome


# Método para buscar el primer en el gen
def __fuzzy_search_single_primer(subseq, seq, mismatches, strand):
    # Based on nt_search from https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/__init__.py
    pattern = ''
    # Se genera el patrón de búsqueda para la expresión regular
    for nt in subseq:
        # Se sustituyen los caracteres IUPAC por sus posibles nucleótidos
        value = IUPACData.ambiguous_dna_values[nt]
        if len(value) == 1:
            pattern += value
        else:
            pattern += '[%s]' % value
    init_pos = -1
    matches = []
    # Se añade el número máximo de mismatches a la expresión regular
    pattern = '({}){{s<={}}}'.format(pattern, mismatches)
    # Se recorren las coincidencias encontradas guardándolas como TargetSequence
    while True:
        init_pos += 1
        # Se elimina la parte ya explorada previamente de la secuencia del gen para mejorar la eficiencia
        s = seq[init_pos:]
        m = re.search(pattern, s)
        if not m:
            break
        current_pos = init_pos
        init_pos += int(m.start(0))
        # Se recoge la secuencia, la posición inicial, final, número de mismatches existentes en la cadena encontrada y la hebra
        matches.append(TargetSequence(
            m.group(), init_pos, m.end() + current_pos, m.fuzzy_counts[0], strand))
    return matches


def reset_genome_batches(genome_batches: list[GenomeBatch]):
    for genome_batch in genome_batches:
        genome_batch.analyzed_genes.clear()
        genome_batch.variants.clear()
        genome_batch.variants_tax_level.clear()
        for genome in genome_batch.genomes:
            genome.amplified_genes.clear()


def trim_genes16s_in_genome(genome: Genome, primer_pair: PrimerPair, output_dir_path: str) -> GenomeBatch:
    forward_matches = __fuzzy_search_single_primer(
        primer_pair.forward_sequence, genome.sequence, MISMATCHES, '+')
    reverse_matches = __fuzzy_search_single_primer(
        primer_pair.reverse_sequence_rc, genome.sequence, MISMATCHES, '+')

    if any(forward_matches) and (reverse_matches):
        for forward_match in forward_matches:
            r_matches = [
                r_match for r_match in reverse_matches if r_match.init_position > forward_match.end_position]
            r_matches_sorted = sorted(r_matches, key=lambda x: x.init_position)
            if any(r_matches_sorted):
                reverse_match = r_matches_sorted[0]
                sequence = genome.sequence[forward_match.init_position:reverse_match.end_position]
                gene = Gene(sequence)
                gene.init_from_params(
                    genome.identifier, forward_match.init_position, reverse_match.end_position, '+')
                if min_size <= len(gene.sequence) <= max_size:
                    genome.amplified_genes.append(gene)

    forward_matches = __fuzzy_search_single_primer(
        primer_pair.forward_sequence, genome.sequence_rev_comp, MISMATCHES, '-')
    reverse_matches = __fuzzy_search_single_primer(
        primer_pair.reverse_sequence_rc, genome.sequence_rev_comp, MISMATCHES, '-')

    if any(forward_matches) and (reverse_matches):
        for forward_match in forward_matches:
            r_matches = [
                r_match for r_match in reverse_matches if r_match.init_position > forward_match.end_position]
            r_matches_sorted = sorted(r_matches, key=lambda x: x.init_position)
            if any(r_matches_sorted):
                reverse_match = r_matches_sorted[0]
                sequence = genome.sequence_rev_comp[forward_match.init_position:reverse_match.end_position]
                gene = Gene(sequence)
                gene.init_from_params(
                    genome.identifier, forward_match.init_position, reverse_match.end_position, '-')
                if min_size <= len(gene.sequence) <= max_size:
                    genome.amplified_genes.append(gene)

    return genome


def trim_genes16s_in_batch(genome_batch: GenomeBatch, primer_pair: PrimerPair, output_dir_path: str) -> GenomeBatch:
    for genome in genome_batch.genomes:
        genome = trim_genes16s_in_genome(genome, primer_pair, output_dir_path)
        genome_batch.analyzed_genes.extend(genome.amplified_genes)
    return genome_batch


def find_variants(genome_batch: GenomeBatch) -> GenomeBatch:
    gene_seqs = [
        gene.sequence for gene in genome_batch.analyzed_genes]
    variant_gene_seqs = set(gene_seqs)
    genome_batch.variants_tax_level = {
        seq: f'v{num}' for num, seq in enumerate(variant_gene_seqs, 1)}

    for gene in genome_batch.analyzed_genes:
        gene.variant_tax_level = genome_batch.variants_tax_level[gene.sequence]
        taxonomy_copy = genome_batch.taxonomy.copy()
        taxonomy_copy['variant'] = gene.variant_tax_level
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


def group_batches_by_tax_level(genome_batches: list[GenomeBatch], taxonomy_level: str) -> dict[str, list[GenomeBatch]]:

    categories_ordered = list(categories_order.keys())
    position = categories_order[taxonomy_level] + 1
    taxonomy_to_group = categories_ordered[:position]

    grouped_batches = {}
    for batch in genome_batches:
        batch_taxonomy = {}
        for tax_level in taxonomy_to_group:
            batch_taxonomy[tax_level] = batch.taxonomy[tax_level]

        batch_taxonomy_str = ';'.join(list(batch_taxonomy.values()))
        if batch_taxonomy_str not in grouped_batches:
            grouped_batches[batch_taxonomy_str] = []
        grouped_batches[batch_taxonomy_str].append(batch)

    return grouped_batches
