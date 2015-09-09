import argparse
import primer3tools

def run():
    parser = argparse.ArgumentParser(
        description = 'Get uniqueness of primers, using outtput of "primer3tools primer3"',
        usage = 'primer3tools primer3 [options] <genomes_file> <primer3_dir> <outprefix>'
    )

    parser.add_argument('--min_product_length', type=int, help='Minimum length of PCR product [%(default)s]', default=50, metavar='INT')
    parser.add_argument('--max_product_length', type=int, help='Maximum length of PCR product [%(default)s]', default=1000, metavar='INT')
    parser.add_argument('genomes_file', help='Input file')
    parser.add_argument('primer3_outdir', help='Primer3 output directory')
    parser.add_argument('outprefix', help='Prefix of output files')

    options = parser.parse_args()

    u = primer3tools.uniqueness.PrimerUniqueness(
        options.genomes_file,
        options.primer3_outdir,
        options.outprefix,
        min_product_length=options.min_product_length,
        max_product_length=options.max_product_length,
    )

    u.run()

