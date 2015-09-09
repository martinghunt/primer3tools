import unittest
import os
from primer3tools import genome_set

modules_dir = os.path.dirname(os.path.abspath(genome_set.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


def write_fake_genomes_file(names, genomes, make_primers, filename):
    assert len(names) == len(genomes) == len(make_primers)
    with open(filename, 'w') as f:
        for i in range(len(genomes)):
            print(names[i], genomes[i], make_primers[i], sep='\t', file=f)



class Genome(unittest.TestCase):
    def test_genome_init(self):
        '''test genome init'''
        fasta_file = os.path.join(data_dir, 'genome_set_dummy_genome.fasta')
        genome = genome_set.Genome(fasta_file, True)
        self.assertEqual(genome.fasta_file, fasta_file)
        self.assertTrue(genome.make_primers)
        genome = genome_set.Genome(fasta_file, False)
        self.assertEqual(genome.fasta_file, fasta_file)
        self.assertFalse(genome.make_primers)

        with self.assertRaises(genome_set.Error):
            genome = genome_set.Genome('thisfiledoesnotexist', True)


class GenomeSet(unittest.TestCase):
    def test_genome_set_init_good_input(self):
        '''test genome_set init with good input'''
        names = ['genome1', 'genome2']
        genomes = ['genome_set_dummy_genome.fasta', 'genome_set_dummy_genome.2.fasta']
        genomes = [os.path.join(data_dir, knightwhosaysni) for knightwhosaysni in genomes]
        config_file = 'tmp.test_genome_set.config'
        write_fake_genomes_file(names, genomes, ['0', '1'], config_file)
        g_set = genome_set.GenomeSet(config_file)
        expected = {
            'genome1': genome_set.Genome(genomes[0], False),
            'genome2': genome_set.Genome(genomes[1], True),
        }
        self.assertEqual(expected, g_set.genomes)
        os.unlink(config_file)


    def test_genome_set_init_bad_input(self):
        '''test genome_set init with bad input'''
        names = ['genome1', 'genome2']
        genomes = ['genome_set_dummy_genome.fasta', 'filedoesnotexist.fa']
        genomes = [os.path.join(data_dir, knightwhosaysni) for knightwhosaysni in genomes]
        config_file = 'tmp.test_genome_set.config'
        write_fake_genomes_file(names, genomes, ['0', '1'], config_file)
        with self.assertRaises(genome_set.Error):
            g_set = genome_set.GenomeSet(config_file)
        os.unlink(config_file)

