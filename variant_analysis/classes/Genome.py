from typing import Optional
from classes.Gene import Gene
import numpy as np
from Bio.Seq import Seq


class Genome:

    header: str
    identifier: str
    taxonomy_identifier: str
    taxonomy: dict[str, str]
    sequence: str
    filename: str
    amplified_genes: list[Gene]
    is_plasmid: bool
    database_name: str
    database_length: int

    def __init__(self, header, identifier, taxonomy_identifier, sequence, filename):
        self.header = header
        self.identifier = identifier
        self.taxonomy_identifier = taxonomy_identifier
        self.taxonomy = {}
        self.sequence = sequence
        self.sequence_rev_comp = str(Seq(sequence).reverse_complement())
        self.filename = filename
        self.amplified_genes = []
        self.is_plasmid = None

    def get_fasta_contents(self):
        return f'{self.header}\n{self.sequence}'

    def get_seq_length(self) -> int:
        return len(self.sequence)


class GenomeBatch:

    taxonomy: dict[str, str]
    taxonomy_str: str
    taxonomy_identifier: str
    genomes: list[Genome]
    genomes_id_str: str
    analyzed_genes: list[Gene]
    variants: dict[str, list[Gene]]
    variants_tax_level: dict[str, str]

    def __init__(self, taxonomy_identifier: str, taxonomy: dict[str, str], genomes: list[Genome]):
        self.taxonomy = taxonomy
        self.taxonomy_str = ';'.join(list(taxonomy.values()))
        self.taxonomy_identifier = taxonomy_identifier
        self.genomes = genomes
        self.genomes_id_str = '-'.join([x.identifier for x in genomes])
        self.analyzed_genes = []
        self.variants = {}
        self.variants_tax_level = {}

    def get_genes_mean_size(self) -> int:
        if self.analyzed_genes:
            mean = np.mean([gene.get_seq_length()
                            for gene in self.analyzed_genes])
        else:
            mean = 0
        return np.around(mean, decimals=3)

    def get_genomes_mean_size(self) -> int:
        mean = np.mean([genome.get_seq_length()
                        for genome in self.genomes])
        return np.around(mean, decimals=3)

    def get_genomes_with_genes_mean_size(self) -> int:
        mean = np.mean([genome.get_seq_length()
                        for genome in self.genomes if genome.amplified_genes])
        return np.around(mean, decimals=3)

    def get_num_analyzed_genes(self) -> int:
        return len(self.analyzed_genes)

    def get_num_variants(self) -> int:
        return len(self.variants)

    def get_num_genes_per_strand(self, strand: str) -> int:
        genes_per_strand = [
            gene for gene in self.analyzed_genes if gene.strand == strand]
        return len(genes_per_strand)

    def get_num_copies_per_variant(self):
        num_copies = {}
        for gene in self.analyzed_genes:
            if gene.updated_taxonomy['variant'] not in num_copies:
                num_copies[gene.updated_taxonomy['variant']] = 0
            num_copies[gene.updated_taxonomy['variant']] += 1
        return num_copies

    def get_main_genome(self) -> Genome:
        main_genome = self.genomes[0]
        for genome in self.genomes:
            if not genome.is_plasmid:
                main_genome = genome
                break
            else:
                main_genome = genome

        return main_genome

    def get_taxonomy_list(self) -> list[str]:
        return list(self.taxonomy.values())

    def get_genomes_with_genes_id_str(self) -> str:
        genomes_with_genes_id_str = '-'.join(
            [x.identifier for x in self.genomes if x.amplified_genes])
        return genomes_with_genes_id_str
