class Gene:

    header: str
    genome_identifier: str
    updated_taxonomy: list[str]
    sequence: str
    init_position: int
    end_position: int
    strand: str
    variant_tax_level: str
    variant_id: str

    def __init__(self, header: str, sequence: str):
        self.header = f'>{header}'
        self.genome_identifier = header.split(";")[0]
        self.updated_taxonomy = []
        self.sequence = sequence
        self.init_position = int(header.split(';')[1].split(',')[0][1:])
        self.end_position = int(header.split(';')[1].split(',')[1][:-2])
        self.strand = header[-2]
        self.variant_tax_level = ''
        self.variant_id = ''

    def get_seq_length(self) -> int:
        return len(self.sequence)
