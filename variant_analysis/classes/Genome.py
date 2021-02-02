from typing import Optional
from classes.Gene import Gene


class Genome:

    header: str
    identifier: str
    taxonomy_identifier: str
    taxonomy: list[str]
    sequence: str
    filename: str
    amplified_genes: list[Gene]
    is_plasmid: Optional[bool]
    database_name: Optional[str]
    database_length: Optional[int]

    def __init__(self, header, identifier, taxonomy_identifier, sequence, filename):
        self.header = header
        self.identifier = identifier
        self.taxonomy_identifier = taxonomy_identifier
        self.taxonomy = []
        self.sequence = sequence
        self.filename = filename
        self.amplified_genes = []
        self.is_plasmid = None
        self.database_name = None
        self.database_length = None

    def get_fasta_contents(self):
        return f'{self.header}\n{self.sequence}'


class GenomeBatch:
    def __init__(self, taxonomy_identifier: str, taxonomy: list[str], genomes: list[Genome]):
        self.taxonomy_identifier = taxonomy_identifier
        self.taxonomy = taxonomy
        self.updated_taxonomy = None
        self.genomes = genomes
        self.genomes_id_str = ', '.join([x.identifier for x in genomes])
        self.genes_to_analyse = []
        self.variant_gene_seqs_ids = {}

    def get_num_genes_to_analyse(self):
        return len(self.genes_to_analyse)

    def get_num_variant_genes(self):
        return len(self.variant_gene_seqs_ids)
