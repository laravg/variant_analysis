from Bio.Seq import Seq


class PrimerPair:
    def __init__(self, identifier, forward_sequence, reverse_sequence):
        self.identifier = identifier
        self.forward_sequence = forward_sequence
        self.reverse_sequence = reverse_sequence
        self.reverse_sequence_rc = str(
            Seq(reverse_sequence).reverse_complement())


# Clase para guardar la información sobre cada secuencia amplificada por un primer
class TargetSequence:
    def __init__(self, sequence, init_position, end_position, mismatches, strand):
        self.sequence = sequence
        self.init_position = init_position
        self.end_position = end_position
        self.mismatches = mismatches
        self.strand = strand
