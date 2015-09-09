import argparse
import primer3tools

def run():
    parser = argparse.ArgumentParser(
        description = 'Run primer3 and bowtie2-build on a batch of genomes',
        usage = 'primer3tools batch [options] <primer3_config> <genomes_file> <outdir>'
    )

    parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1)
    parser.add_argument('primer3_config', help='Primer3 config file')
    parser.add_argument('genomes_file', help='File of genomes information')
    parser.add_argument('outdir', help='Primer3 output directory')
    options = parser.parse_args()

    batch = primer3tools.primer3_batch.Primer3Batch(
        options.primer3_config,
        options.genomes_file,
        options.outdir,
        threads=options.threads
    )
    batch.run()
