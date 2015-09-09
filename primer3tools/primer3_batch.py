import os
import multiprocessing
import primer3tools


def _make_directory(d):
    if os.path.exists(d):
        return

    try:
        os.mkdir(d)
    except:
        raise Error('Error mkdir ' + d)


def _run_analysis(genome_name, genome, primer3_config, outprefix):
    primers_fasta = outprefix + '.primers.fasta.gz'
    primer3_outfile = outprefix + '.primer3_core.out.gz'

    if genome.make_primers and ((not os.path.exists(primers_fasta)) or (not os.path.exists(primer3_outfile))):
        p3 = primer3tools.primer3.Primer3(genome.fasta_file, primer3_config, genome_name)
        p3.run(outprefix)

    bowtie2_index_dir = outprefix + '.bowtie2_index'
    bowtie2_index = os.path.join(bowtie2_index_dir, 'index')

    if not primer3tools.mapping.is_bowtie2_indexed(bowtie2_index):
        _make_directory(bowtie2_index_dir)
        primer3tools.mapping.bowtie2_index(genome.fasta_file, outprefix=bowtie2_index)


# throws a pickle error without this wrapper...
def _run_analysis_wrapper(y):
    _run_analysis(*y)
    return y


class Error (Exception): pass


class Primer3Batch:
    def __init__(self, primer3_config, genomes_file, primer3_outdir, threads=1):
        self.primer3_config = os.path.abspath(primer3_config)
        self.genomes_file = os.path.abspath(genomes_file)
        self.primer3_outdir = os.path.abspath(primer3_outdir)
        self.threads = threads


    def run(self):
        genomes = primer3tools.genome_set.GenomeSet(self.genomes_file)
        _make_directory(self.primer3_outdir)
        pool = multiprocessing.Pool(self.threads)

        x = [(name, genomes[name], self.primer3_config, os.path.join(self.primer3_outdir, name)) for name in genomes]
        pool.map(_run_analysis_wrapper, x)

