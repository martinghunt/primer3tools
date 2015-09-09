import os


class Error (Exception): pass


class Genome:
    def __init__(self, fasta_file, make_primers):
        self.fasta_file = os.path.abspath(fasta_file)
        self.make_primers = make_primers
        if not os.path.exists(self.fasta_file):
            raise Error('Genome fasta file not found: ' + self.fasta_file)


    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__


    def __str__(self):
        return self.fasta_file + '\t' + str(self.make_primers)


class GenomeSet:
    def __init__(self, infile):
        self.infile = os.path.abspath(infile)
        self.genomes = self._parse_input_file(self.infile)


    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__ and self.genomes == other.genomes


    def _parse_input_file(self, infile):
        genomes = {}


        with open(infile) as f:
            for line in f:
                genome_name, filename, make_primers = line.rstrip().split('\t')
                if genome_name in genomes:
                   raise Error('Non-unique genome name: ' + genome_name)
                if make_primers not in ['0', '1']:
                   raise Error('Final column of input file must be 0 or 1. Bad line:\n' + line)
                genomes[genome_name] = Genome(filename, make_primers == '1')

        return genomes


    def __getitem__(self, key):
        return self.genomes[key]


    def __iter__(self):
        for key in self.genomes:
            yield key
