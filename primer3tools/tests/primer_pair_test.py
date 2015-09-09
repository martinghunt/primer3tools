import unittest
import os
from primer3tools import primer_pair
from pyfastaq import sequences

modules_dir = os.path.dirname(os.path.abspath(primer_pair.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestPrimer3(unittest.TestCase):
    def test_primer_pair_init_fails(self):
        '''test PrimerPair __init__ fails'''
        data_dict = {}
        new_keys_and_values = [
            ('SEQUENCE_ID', 'name'),
            ('PRIMER_LEFT_0', '42,4'),
            ('PRIMER_RIGHT_0', '100,5'),
            ('PRIMER_LEFT_0_SEQUENCE', 'ACGT'),
            ('PRIMER_RIGHT_0_SEQUENCE', 'ACGTA'),
        ]

        for key, value in new_keys_and_values:
            with self.assertRaises(primer_pair.Error):
                pair = primer_pair.PrimerPair(data_dict, 0, 'genome_name')
            data_dict[key] = value

        pair = primer_pair.PrimerPair(data_dict, 0, 'genome_name')
        with self.assertRaises(primer_pair.Error):
            pair = primer_pair.PrimerPair(data_dict, 1, 'genome_name')


    def test_primer_pair_init_passes(self):
        '''test PrimerPair __init__ passes'''
        input_dict = {
            'SEQUENCE_ID': 'seq_name',
            'PRIMER_LEFT_0': '42,4',
            'PRIMER_RIGHT_0': '100,4',
            'PRIMER_LEFT_0_SEQUENCE': 'AAAA',
            'PRIMER_RIGHT_0_SEQUENCE': 'AGGT',
        }


        left_fasta = sequences.Fasta('genome_name__seq_name__0__42__97/1', 'AAAA')
        right_fasta = sequences.Fasta('genome_name__seq_name__0__42__97/2', 'AGGT')
        pair = primer_pair.PrimerPair(input_dict, 0, 'genome_name')
        self.assertEqual(pair.left_fasta, left_fasta)
        self.assertEqual(pair.right_fasta, right_fasta)

