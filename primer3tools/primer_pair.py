from pyfastaq import sequences

class Error (Exception): pass


class PrimerPair:
    def __init__(self, data_dict, index, genome_name):
        self.index = index
        self.genome_name = genome_name
        self.sequence_id = None
        self.left_start = None
        self.right_start = None
        self.left_fasta = None
        self.right_fasta = None
        self.hits_to_genomes = {}

        if data_dict is not None and self.index is not None:
            try:
                self.sequence_id = data_dict.get('SEQUENCE_ID', None)

                # primer3 reports the end position of the sequence, and it is reverse complemented,
                # so the primer pairs look like proper read pairs
                left_sequence = data_dict.get('PRIMER_LEFT_' + str(self.index) + '_SEQUENCE', None)
                right_sequence = data_dict.get('PRIMER_RIGHT_' + str(self.index) + '_SEQUENCE', None)
                self.left_start = int(data_dict['PRIMER_LEFT_' + str(self.index)].split(',')[0])
                right_end = int(data_dict['PRIMER_RIGHT_' + str(self.index)].split(',')[0])
                self.right_start = right_end - len(right_sequence) + 1

                name_prefix = '__'.join([self.genome_name, self.sequence_id, str(self.index), str(self.left_start), str(self.right_start)])
                self.left_fasta = sequences.Fasta(name_prefix + '/1', left_sequence)
                self.right_fasta = sequences.Fasta(name_prefix + '/2', right_sequence)
            except:
                raise Error('Error making PrimerPair from data_dict:\n' + str(data_dict))


        if not self.has_all_info():
            raise Error('Error making PrimerPair from data_dict:\n' + str(data_dict))


    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__


    def has_all_info(self):
        return None not in [
            self.sequence_id,
            self.left_start,
            self.right_start,
            self.left_fasta,
            self.right_fasta
        ]


    def __str__(self):
        return '\t'.join([
            self.sequence_id,
            str(self.index),
            str(self.left_start),
            str(self.right_start),
            self.left_fasta.seq,
            self.right_fasta.seq
        ])

