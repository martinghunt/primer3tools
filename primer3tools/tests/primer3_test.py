import unittest
import pyfastaq
import shutil
import filecmp
from unittest.mock import patch
import os
from primer3tools import primer3, primer_pair

modules_dir = os.path.dirname(os.path.abspath(primer3.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestPrimer3(unittest.TestCase):
    def setUp(self):
        dummy_config = os.path.join(data_dir, 'primer3_test_dummy.config')
        dummy_fasta = os.path.join(data_dir, 'primer3_test_dummy.fa')
        self.p3 = primer3.Primer3(dummy_fasta, dummy_config, 'genome_name')


    def test_split_primer3_output_line(self):
        '''test _split_primer3_output_line'''
        with self.assertRaises(primer3.Error):
            self.p3._split_primer3_output_line('does not have an equals sign')

        tests = [
            ('key=value', 'key', 'value'),
            ('key=value\n', 'key', 'value'),
            ('key=value=foo', 'key', 'value=foo'),
            ('key=value=foo\n', 'key', 'value=foo'),
        ]

        for line, key, value in tests:
            self.assertEqual(self.p3._split_primer3_output_line(line), (key, value))


    def test_get_next_primer3_sequence_results(self):
        '''test _get_next_primer3_sequence_results'''
        infile = os.path.join(data_dir, 'primer3_test_get_next_primer3_sequence_results.infile')
        expected_dicts = [
            {'ONE1': 'one1', 'ONE2': 'one2'},
            {'TWO1': 'two1'},
            None
        ]

        with open(infile) as f:
            for expected in expected_dicts:
                got = self.p3._get_next_primer3_sequence_results(f)
                self.assertEqual(expected, got)


    def test_get_next_primer3_sequence_results_fails(self):
        '''test _get_next_primer3_sequence_results on bad input'''
        infile = os.path.join(data_dir, 'primer3_test_get_next_primer3_sequence_results.infile_bad')
        with open(infile) as f:
            self.p3._get_next_primer3_sequence_results(f)
            with self.assertRaises(primer3.Error):
                self.p3._get_next_primer3_sequence_results(f)


    def test_primer3_sequence_results_to_list(self):
        '''test _primer3_sequence_results_to_list'''
        pair1_dict = {
            'SEQUENCE_ID': 'name_of_seq',
            'PRIMER_LEFT_0': '42,4',
            'PRIMER_RIGHT_0': '100,4',
            'PRIMER_LEFT_0_SEQUENCE': 'AAAA',
            'PRIMER_RIGHT_0_SEQUENCE': 'CCCC',
        }

        pair2_dict = {
            'SEQUENCE_ID': 'name_of_seq',
            'PRIMER_LEFT_1': '42,3',
            'PRIMER_RIGHT_1': '100,3',
            'PRIMER_LEFT_1_SEQUENCE': 'GGG',
            'PRIMER_RIGHT_1_SEQUENCE': 'TTT',
        }

        input_dict = {
            'PRIMER_LEFT_NUM_RETURNED': 2,
            'PRIMER_RIGHT_NUM_RETURNED': 2,
            'PRIMER_PAIR_NUM_RETURNED': 2,
            'SEQUENCE_ID': 'name_of_seq',
        }

        input_dict.update(pair1_dict)

        with self.assertRaises(primer_pair.Error):
            self.p3._primer3_sequence_results_to_list(input_dict)

        input_dict.update(pair2_dict)
        got = self.p3._primer3_sequence_results_to_list(input_dict)

        expected = [
            primer_pair.PrimerPair(pair1_dict, 0, 'genome_name'),
            primer_pair.PrimerPair(pair2_dict, 1, 'genome_name'),
        ]

        self.assertEqual(expected, got)


    def test_has_no_primers(self):
        '''test _has_no_primers'''
        data_dicts = [
            {'PRIMER_PAIR_NUM_RETURNED': '0'},
            {'PRIMER_ERROR': 'SEQUENCE_INCLUDED_REGION length < min PRIMER_PRODUCT_SIZE_RANGE'}
        ]

        for data_dict in data_dicts:
            self.assertEqual(True, self.p3._primer3_found_no_primers(data_dict))


    def test_primer3_sequence_results_to_list_fails(self):
        '''test _primer3_sequence_results_to_list on bad input'''
        inputs = [
            {},
            {'PRIMER_LEFT_NUM_RETURNED': 1, 'PRIMER_RIGHT_NUM_RETURNED': 1},
            {'PRIMER_LEFT_NUM_RETURNED': 1, 'PRIMER_RIGHT_NUM_RETURNED': 1, 'PRIMER_PAIR_NUM_RETURNED': 2},
            {'PRIMER_RIGHT_NUM_RETURNED': 1, 'PRIMER_PAIR_NUM_RETURNED': 1},
            {'PRIMER_LEFT_NUM_RETURNED': 1, 'PRIMER_PAIR_NUM_RETURNED': 1},
        ]

        for results in inputs:
            with self.assertRaises(primer3.Error):
                self.p3._primer3_sequence_results_to_list(results)

        makes_primer_pair_error = {'PRIMER_LEFT_NUM_RETURNED': 1, 'PRIMER_RIGHT_NUM_RETURNED': 1, 'PRIMER_PAIR_NUM_RETURNED': 1}

        with self.assertRaises(primer_pair.Error):
            self.p3._primer3_sequence_results_to_list(makes_primer_pair_error)


    def test_load_primer_pairs(self):
        '''test _load_primer_pairs'''
        infile = os.path.join(data_dir, 'primer3_test_load_primer_pairs.infile')
        got = self.p3._load_primer_pairs(infile)

        seq1_dict = {
            'SEQUENCE_ID': 'seq1',
            'PRIMER_LEFT_NUM_RETURNED': '1',
            'PRIMER_RIGHT_NUM_RETURNED': '1',
            'PRIMER_PAIR_NUM_RETURNED': '1',
            'PRIMER_LEFT_0': '42,4',
            'PRIMER_RIGHT_0': '100,3',
            'PRIMER_LEFT_0_SEQUENCE': 'AAAA',
            'PRIMER_RIGHT_0_SEQUENCE': 'CCC',
        }

        seq2_dict = {
            'SEQUENCE_ID': 'seq2',
            'PRIMER_LEFT_NUM_RETURNED': '2',
            'PRIMER_RIGHT_NUM_RETURNED': '2',
            'PRIMER_PAIR_NUM_RETURNED': '2',
            'PRIMER_LEFT_0': '142,2',
            'PRIMER_RIGHT_0': '200,1',
            'PRIMER_LEFT_0_SEQUENCE': 'AA',
            'PRIMER_RIGHT_0_SEQUENCE': 'C',
            'PRIMER_LEFT_1': '242,4',
            'PRIMER_RIGHT_1': '300,4',
            'PRIMER_LEFT_1_SEQUENCE': 'GGGG',
            'PRIMER_RIGHT_1_SEQUENCE': 'TTTT',
        }

        expected = {
            'seq1': [primer_pair.PrimerPair(seq1_dict, 0, 'genome_name')],
            'seq2': [
                primer_pair.PrimerPair(seq2_dict, 0, 'genome_name'),
                primer_pair.PrimerPair(seq2_dict, 1, 'genome_name'),
            ]
        }

        self.assertEqual(expected, got)


    def test_write_primers_fasta(self):
        '''test _write_primers_fasta'''
        seq1_dict = {
            'SEQUENCE_ID': 'seq1',
            'PRIMER_LEFT_NUM_RETURNED': '1',
            'PRIMER_RIGHT_NUM_RETURNED': '1',
            'PRIMER_PAIR_NUM_RETURNED': '1',
            'PRIMER_LEFT_0': '42,4',
            'PRIMER_RIGHT_0': '100,3',
            'PRIMER_LEFT_0_SEQUENCE': 'AAAA',
            'PRIMER_RIGHT_0_SEQUENCE': 'CCC',
        }

        seq2_dict = {
            'SEQUENCE_ID': 'seq2',
            'PRIMER_LEFT_NUM_RETURNED': '2',
            'PRIMER_RIGHT_NUM_RETURNED': '2',
            'PRIMER_PAIR_NUM_RETURNED': '2',
            'PRIMER_LEFT_0': '142,2',
            'PRIMER_RIGHT_0': '200,1',
            'PRIMER_LEFT_0_SEQUENCE': 'AA',
            'PRIMER_RIGHT_0_SEQUENCE': 'C',
            'PRIMER_LEFT_1': '242,4',
            'PRIMER_RIGHT_1': '300,4',
            'PRIMER_LEFT_1_SEQUENCE': 'GGGG',
            'PRIMER_RIGHT_1_SEQUENCE': 'TTTT',
        }
        self.p3.primer_pairs = {
            'seq1': [primer_pair.PrimerPair(seq1_dict, 0, 'genome_name1')],
            'seq2': [
                primer_pair.PrimerPair(seq2_dict, 0, 'genome_name2'),
                primer_pair.PrimerPair(seq2_dict, 1, 'genome_name2'),
            ]
        }
        outfile = 'tmp.primer_pairs.fa'
        self.p3._write_primers_fasta(outfile)
        expected = os.path.join(data_dir, 'primer3_test_write_primers_fasta.out.fa')
        self.assertTrue(filecmp.cmp(outfile, expected, shallow=False))
        os.unlink(outfile)


    @patch('primer3tools.primer3.Primer3._run_primer3_core')
    def test_run(self, primer3_mock):
        '''test run'''
        original_primer3_outfile = os.path.join(data_dir, 'primer3_test_run.primer3_core.out.gz')
        primer3_outprefix = 'tmp.primer3_test_run'
        primer3_outfile = primer3_outprefix + '.primer3_core.out.gz'
        shutil.copy(original_primer3_outfile, primer3_outfile)
        fasta_outfile = primer3_outprefix + '.primers.fasta.gz'
        expected_outfile = os.path.join(data_dir, 'primer3_test_run.expected.out.fa.gz')
        expected_seqs = {}
        pyfastaq.tasks.file_to_dict(expected_outfile, expected_seqs)
        self.p3.run(primer3_outprefix)
        got_seqs = {}
        pyfastaq.tasks.file_to_dict(fasta_outfile, got_seqs)
        self.assertEqual(got_seqs, expected_seqs)
        os.unlink(fasta_outfile)
        os.unlink(primer3_outfile)

