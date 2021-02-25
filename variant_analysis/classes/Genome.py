from typing import Optional
from classes.Gene import Gene


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
        self.filename = filename
        self.amplified_genes = []

    def get_fasta_contents(self):
        return f'{self.header}\n{self.sequence}'


class GenomeBatch:

    taxonomy: str
    taxonomy_dict: dict[str, str]
    taxonomy_ids: list[str]
    updated_taxonomy: dict[str, str]
    genomes: list[Genome]
    genomes_id_str: str
    analyzed_genes: list[Gene]
    variants: dict[str, list[Gene]]
    variants_tax_level: dict[str, str]

    def __init__(self, taxonomy: str, taxonomy_dict: dict[str, str], taxonomy_ids: list[list[str]], genomes: list[Genome]):
        self.taxonomy = taxonomy
        self.taxonomy_dict = taxonomy_dict
        self.taxonomy_ids = taxonomy_ids
        self.updated_taxonomy = []
        self.genomes = genomes
        self.genomes_id_str = '_'.join([x.identifier for x in genomes])
        self.analyzed_genes = []
        self.variants = {}
        self.variants_tax_level = {}

    def get_num_analyzed_genes(self) -> int:
        return len(self.analyzed_genes)

    def get_num_variants(self) -> int:
        return len(self.variants)

    def get_main_genome(self) -> Genome:
        main_genome = self.genomes[0]
        for genome in self.genomes:
            if not genome.is_plasmid:
                main_genome = genome
                break
            else:
                main_genome = genome

        return main_genome
