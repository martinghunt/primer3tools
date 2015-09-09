import os
from primer3tools import common


class Error (Exception): pass


bowtie2_index_extensions = [x + '.bt2' for x in ['1', '2', '3', '4', 'rev.1', 'rev.2']]


def is_bowtie2_indexed(infile):
    for ext in bowtie2_index_extensions:
        if not os.path.exists(infile + '.' + ext):
            return False

    return True


def bowtie2_index(infile, outprefix=None):
    if outprefix is None:
        outprefix = infile

    if not is_bowtie2_indexed(infile):
        common.syscall('bowtie2-build ' + infile + ' ' + outprefix)


def run_bowtie2(reads, reference, outfile, threads=1):
    assert is_bowtie2_indexed(reference)
    cmd = ' '.join([
        'bowtie2',
        '-x', reference,
        '-U', reads,
        '-f',  # reads are in fasta format
        '--end-to-end --very-fast',
        '--threads', str(threads),
        '--all', # report all alignments
        '--reorder', # force SAM output order to match order of input reads
        '-S', outfile # output in SAM format
    ])

    common.syscall(cmd)




