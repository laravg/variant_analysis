class Gene:
    def __init__(self, header, sequence):
        self.header = f'>{header}'
        self.genome_identifier = header.split(";")[0]
        self.updated_taxonomy = None
        self.sequence = sequence
        self.init_position = header.split(';')[1].split(',')[0].strip('[')
        self.end_position = header.split(';')[1].split(',')[1].strip(']')
        self.strand = header[-2]
        self.variant_id = None
