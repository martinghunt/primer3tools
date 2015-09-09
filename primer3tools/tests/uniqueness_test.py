import unittest
import pysam
import filecmp
import os
from primer3tools import uniqueness, genome_set

modules_dir = os.path.dirname(os.path.abspath(uniqueness.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


def write_fake_genomes_file(names, genomes, make_primers, filename):
    assert len(names) == len(genomes) == len(make_primers)
    with open(filename, 'w') as f:
        for i in range(len(genomes)):
            print(names[i], genomes[i], make_primers[i], sep='\t', file=f)


class TestPrimer3(unittest.TestCase):
    def test_cat_primer_fastas(self):
        '''test _cat_primer_fastas'''
        primer3_dir = os.path.join(data_dir, 'uniqueness_test_cat_primer_fastas.primer3_dir')
        names = ['genome1', 'genome2', 'genome3']
        genome_files = [
            os.path.join(data_dir, 'uniqueness_test_cat_primer_fastas.genome1.fa'),
            os.path.join(data_dir, 'uniqueness_test_cat_primer_fastas.genome2.fa'),
            os.path.join(data_dir, 'uniqueness_test_cat_primer_fastas.genome3.fa')
        ]
        make_primers = [1, 1, 0]
        genomes_file = 'tmp.test.uniqueness_cat_primer_fastas.genomes_file'
        write_fake_genomes_file(names, genome_files, make_primers, genomes_file)
        uniq = uniqueness.PrimerUniqueness(genomes_file, primer3_dir, 'outprefix')
        genomes = genome_set.GenomeSet(genomes_file)
        catted_fasta = 'tmp.test.uniqueness_cat_primer_fastas.out.fa'
        uniq._cat_primer_fastas(genomes, catted_fasta)
        expected = os.path.join(data_dir, 'uniqueness_test_cat_primer_fastas.out.fa')
        self.assertTrue(filecmp.cmp(expected, catted_fasta, shallow=False))
        os.unlink(genomes_file)
        os.unlink(catted_fasta)


    def test_cat_all_genomes(self):
        '''test _cat_all_genomes'''
        names = ['genome1', 'genome2', 'genome3']
        genome_files = [
            os.path.join(data_dir, 'uniqueness_test_cat_all_genomes.genome1.fa'),
            os.path.join(data_dir, 'uniqueness_test_cat_all_genomes.genome2.fa'),
            os.path.join(data_dir, 'uniqueness_test_cat_all_genomes.genome3.fa')
        ]
        make_primers = [1, 1, 0]
        genomes_file = 'tmp.test.uniqueness_cat_all_genomes.genomes_file'
        write_fake_genomes_file(names, genome_files, make_primers, genomes_file)
        genomes = genome_set.GenomeSet(genomes_file)
        catted_fasta = 'tmp.test.uniqueness_cat_all_genomes.out.fa'
        uniqueness.PrimerUniqueness._cat_all_genomes(genomes, catted_fasta)
        expected = os.path.join(data_dir, 'uniqueness_test_cat_all_genomes.out.fa')
        self.assertTrue(filecmp.cmp(expected, catted_fasta, shallow=False))
        os.unlink(genomes_file)
        os.unlink(catted_fasta)


    def test_is_perfect_hit(self):
        '''test _is_perfect_hit'''
        sam_file = os.path.join(data_dir, 'uniqueness_test_is_perfect_hit.sam')
        sam_reader = pysam.Samfile(sam_file, "r")
        sam_records = [record for record in sam_reader.fetch(until_eof=True)]
        self.assertFalse(uniqueness.PrimerUniqueness._is_perfect_hit(sam_records[0]))
        self.assertTrue(uniqueness.PrimerUniqueness._is_perfect_hit(sam_records[1]))
        self.assertFalse(uniqueness.PrimerUniqueness._is_perfect_hit(sam_records[2]))
        self.assertFalse(uniqueness.PrimerUniqueness._is_perfect_hit(sam_records[4]))
        with self.assertRaises(uniqueness.Error):
            uniqueness.PrimerUniqueness._is_perfect_hit(sam_records[3])


    def test_sam_to_read_sequence(self):
        '''test _sam_to_read_sequence'''
        expected = ['AATGTGTATAGCGT', 'ATATTCGGCCACGG']
        sam_file = os.path.join(data_dir, 'uniqueness_test_sam_to_read_sequence.sam')
        sam_reader = pysam.Samfile(sam_file, "r")
        sam_records = [record for record in sam_reader.fetch(until_eof=True)]
        assert len(sam_records) == len(expected)
        for i in range(len(sam_records)):
            self.assertEqual(expected[i], uniqueness.PrimerUniqueness._sam_to_read_sequence(sam_records[i]))


    def test_filter_pairs_dict(self):
        '''test _filter_pairs_dict'''
        test_dict = {
            'primer1': {'left': {}, 'right': {'contig1': ['hit1']}},
            'primer2': {'right': {}, 'left': {'contig1': ['hit1']}},
            'primer3': {'left': {'contig2': ['hit1']}, 'right': {'contig1': ['hit1']}},
        }
        expected = {
            'primer3': {'left': {'contig2': ['hit1']}, 'right': {'contig1': ['hit1']}},
        }

        uniqueness.PrimerUniqueness._filter_pairs_dict(test_dict)
        self.assertEqual(test_dict, expected)


    def test_parse_sam(self):
        '''test _parse_sam'''
        sam_file = os.path.join(data_dir, 'uniqueness_test_parse_sam.sam')
        got_dict = uniqueness.PrimerUniqueness._parse_sam(sam_file)
        expected_dict = {
            'primer2': {
                'left': {'contig1': [(0, 13, False, 'AGTAATTAATAAC')]},
                'right': {'contig2': [(31, 14, True, 'TCGCTCCAGGTACG')]}
            }
        }
        self.assertEqual(got_dict, expected_dict)


    def test_is_good_primer_pair(self):
        '''test _is_good_primer_pair'''
        tests = [
            [False, (0, 10, False, 'ACGTACGTAC'), (30, 10, False, 'ACGTACGTAC')], # wrong orientation
            [False, (0, 10, True, 'ACGTACGTAC'), (30, 10, True, 'ACGTACGTAC')], # wrong orientation
            [False, (0, 10, False, 'ACGTACGTAC'), (30, 10, True, 'ACGTACGTAC')], # wrong orientation
            [False, (0, 10, False, 'ACGTACGTAC'), (39, 10, True, 'ACGTACGTAC')],  # too short
            [False, (0, 10, False, 'ACGTACGTAC'), (500, 10, True, 'ACGTACGTAC')],  # too long
            [False, (39, 10, True, 'ACGTACGTAC'), (0, 10, False, 'ACGTACGTAC')],  # too short
            [False, (30, 10, False, 'ACGTACGTAC'), (0, 10, True, 'ACGTACGTAC')],  # wrong order
            [True, (0, 10, False, 'ACGTACGTAC'), (40, 10, True, 'ACGTACGTAC')],
            [True, (40, 10, True, 'ACGTACGTAC'), (0, 10, False, 'ACGTACGTAC')],
        ]

        for expected, left_primer, right_primer in tests:
            u = uniqueness.PrimerUniqueness('genomes_file', 'primer3_outdir', 'outprefix', min_product_length=50, max_product_length=200)
            self.assertEqual(expected, u._is_good_primer_pair(left_primer, right_primer))


    def test_good_primer_pairs_from_lists(self):
        '''test _good_primer_pairs_from_lists'''
        p1 = (0, 10, False, 'ACGTACGTAC')
        p2 = (30, 10, False, 'ACGTACGTAC')
        p3 = (40, 10, True, 'ACGTACGTAC')
        p4 = (35, 10, True, 'ACGTACGTAC')
        p5 = (1000, 10, False, 'ACGTACGTAC')
        p6 = (1040, 10, True, 'ACGTACGTAC')

        list1 = [p1, p5]
        list2 = [p2, p3, p4, p6]
        expected = [(p1, p3), (p5, p6)]

        u = uniqueness.PrimerUniqueness('genomes_file', 'primer3_outdir', 'outprefix', min_product_length=50, max_product_length=200)
        got = u._good_primer_pairs_from_lists(list1, list2)
        self.assertEqual(expected, got)


    def test_all_primer_matches(self):
        '''test _all_primer_matches'''
        test_dict = {
            'left': {
                'contig1': [(0, 10, False, 'ACGTACGTAC'), (40, 10, True, 'ACGTACGTAC')],
                'contig2': [(0, 10, False, 'ACGTACGTAC'), (1000, 10, False, 'ACGTACGTAC')],
            },
            'right': {
                'contig1': [(40, 10, True, 'ACGTACGTAC')],
                'contig2': [(1040, 10, True, 'ACGTACGTAC')],
            },
        }

        u = uniqueness.PrimerUniqueness('genomes_file', 'primer3_outdir', 'outprefix', min_product_length=50, max_product_length=200)
        got = u._all_primer_matches(test_dict)
        expected = {
            'contig1': [((0, 10, False, 'ACGTACGTAC'), (40, 10, True, 'ACGTACGTAC'))],
            'contig2': [((1000, 10, False, 'ACGTACGTAC'), (1040, 10, True, 'ACGTACGTAC'))]
        }
        self.assertEqual(expected, got)


    def test_update_primer_hits(self):
        '''test _update_primer_hits'''
        u = uniqueness.PrimerUniqueness('genomes_file', 'primer3_outdir', 'outprefix', min_product_length=50, max_product_length=200)
        primer_hits = {}
        u._update_primer_hits(primer_hits, {}, 'genome1')
        self.assertEqual({}, primer_hits)

        pairs_dict = {
            'primer1': {
                'left': {'contig1': [(0, 13, False, 'AGTAATTAATAAC')]},
                'right': {'contig1': [(100, 14, True, 'TCGCTCCAGGTACG')]}
            }
        }
        u._update_primer_hits(primer_hits, pairs_dict, 'genome2')

        expected = {
            'primer1': {
                'genome2': {
                    'contig1': [((0, 13, False, 'AGTAATTAATAAC'), (100, 14, True, 'TCGCTCCAGGTACG'))]
                }
            }
        }
        self.assertEqual(expected, primer_hits)


    def test_write_all_output_files(self):
        '''test _write_all_output_files'''
        primer3_dir = os.path.join(data_dir, 'uniqueness_test_cat_primer_fastas.primer3_dir')
        names = ['genome1', 'genome2', 'genome3']
        genome_files = [
            os.path.join(data_dir, 'uniqueness_test_cat_primer_fastas.genome1.fa'),
            os.path.join(data_dir, 'uniqueness_test_cat_primer_fastas.genome2.fa'),
            os.path.join(data_dir, 'uniqueness_test_cat_primer_fastas.genome3.fa')
        ]
        make_primers = [1, 1, 0]
        genomes_file = 'tmp.test.uniqueness_cat_primer_fastas.genomes_file'
        write_fake_genomes_file(names, genome_files, make_primers, genomes_file)
        uniq = uniqueness.PrimerUniqueness(genomes_file, primer3_dir, 'tmp.test_write_all_output_files')
        genomes = genome_set.GenomeSet(genomes_file)

        primer_hits = {
            'primer1': {
                'genome2': {
                    'contig1': [((0, 13, False, 'AGTAATTAATAAC'), (100, 14, True, 'TCGCTCCAGGTACG'))]
                }
            },
            'primer2': {
                'genome1': {
                    'contig1': [((0, 13, False, 'AGTAATTAATAAC'), (100, 14, True, 'TCGCTCCAGGTACG'))]
                },
                'genome2': {
                    'contig2': [((10, 23, False, 'AGTAATTAATAAC'), (110, 24, True, 'TCGCTCCAGGTACG'))]
                }
            }
        }
        uniq._write_all_output_files(primer_hits, genomes)
        for suffix in ('all_primers.hits.tsv', 'genome_uniqueness.tsv', 'unique_primers.tsv'):
            got = 'tmp.test_write_all_output_files.' + suffix
            expected = os.path.join(data_dir, 'uniqueness_test_write_all_output_files.expected.' + suffix)
            self.assertTrue(filecmp.cmp(got, expected, shallow=False))
            os.unlink(got)
        os.unlink(genomes_file)


    def test_is_reverse_to_string(self):
        '''test _is_reverse_to_string'''
        self.assertEqual(uniqueness.PrimerUniqueness._is_reverse_to_string(True), '-')
        self.assertEqual(uniqueness.PrimerUniqueness._is_reverse_to_string(False), '+')

