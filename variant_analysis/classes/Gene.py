class Gene:

    header: str
    genome_identifier: str
    updated_taxonomy: dict[str, str]
    sequence: str
    init_position: int
    end_position: int
    strand: str
    variant_tax_level: str
    variant_id: str

    def __init__(self, sequence: str):
        self.updated_taxonomy = []
        self.sequence = sequence
        self.variant_tax_level = ''
        self.variant_id = ''

    def init_from_header(self, header):
        self.header = f'>{header}'
        self.genome_identifier = header.split(";")[0]
        self.init_position = int(header.split(';')[1].split(',')[0][1:])
        self.end_position = int(header.split(';')[1].split(',')[1][:-2])
        self.strand = header[-2]

    def init_from_params(self, genome_identifier, init_position, end_position, strand):
        self.genome_identifier = genome_identifier
        self.init_position = init_position
        self.end_position = end_position
        self.strand = strand

    def get_seq_length(self) -> int:
        return len(self.sequence)

    def get_updated_taxonomy_list(self) -> list[str]:
        return list(self.updated_taxonomy.values())

    def get_updated_taxonomy_str(self) -> str:
        return ";".join(self.get_updated_taxonomy_list())
