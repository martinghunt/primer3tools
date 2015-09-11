import os
import tempfile
import shutil
import pyfastaq
from primer3tools import common, primer_pair


class Error (Exception): pass

class Primer3:
    def __init__(self, fasta_file, config_file, genome_name, primer3_command='primer3_core'):
        self.input_fasta = os.path.abspath(fasta_file)
        self.config_file = os.path.abspath(config_file)
        self.genome_name = genome_name
        self.primer3_command = primer3_command


        for filename in [self.input_fasta, self.config_file]:
            if not os.path.exists(filename):
                raise Error('File not found: "' + filename + '". Cannot continue')


        if shutil.which(self.primer3_command) is None:
            raise Error('Error: primer3 command not found: ' + self.primer3_command)


    def _run_primer3_core(self, fasta_file, config_file, outfile):
        tmpdir = tempfile.mkdtemp(prefix='tmp.run_primer3_core.', dir=os.getcwd())
        boulder_file = os.path.join(tmpdir, 'in.boulder')
        pyfastaq.tasks.to_boulderio(self.input_fasta, boulder_file)
        cmd = ' '.join([
            self.primer3_command,
            '-p3_settings_file=' + self.config_file,
            '<',
            boulder_file,
            '| gzip -9 -c >',
            outfile,
        ])

        common.syscall(cmd)
        shutil.rmtree(tmpdir)


    def _split_primer3_output_line(self, line):
        if '=' not in line:
            raise Error('Error parsing primer3_core output, no equals in this line:' + line)

        return tuple(line.rstrip().split('=', maxsplit=1))


    def _get_next_primer3_sequence_results(self, filehandle):
        results = {}
        line = filehandle.readline()
        if line == '':
            return None

        while line != '=\n':
            if line == '':
                raise Error('Error reading primer3_core output file. Got to end file before a line of just "="')

            key, value = self._split_primer3_output_line(line)
            assert key not in results
            results[key] = value
            line = filehandle.readline()

        return results


    def _primer3_sequence_results_to_list(self, results):
        try:
            number_of_primer_pairs_left = int(results['PRIMER_LEFT_NUM_RETURNED'])
            number_of_primer_pairs_right = int(results['PRIMER_RIGHT_NUM_RETURNED'])
            number_of_primer_pairs = int(results['PRIMER_PAIR_NUM_RETURNED'])
            if len({number_of_primer_pairs_left, number_of_primer_pairs_right, number_of_primer_pairs}) != 1:
                raise('Mismatching number of pairs from primer3_core output file')
        except:
            raise Error('Error getting number of primer pairs from primer3_core output file')

        primer_pairs = []

        for i in range(number_of_primer_pairs):
            primer_pairs.append(primer_pair.PrimerPair(results, i, self.genome_name))

        return primer_pairs


    def _primer3_found_no_primers(self, data_dict):
        return (
          ('PRIMER_ERROR' in data_dict and data_dict['PRIMER_ERROR'] == 'SEQUENCE_INCLUDED_REGION length < min PRIMER_PRODUCT_SIZE_RANGE') or
          ('PRIMER_PAIR_NUM_RETURNED' in data_dict and data_dict['PRIMER_PAIR_NUM_RETURNED'] == '0')
        )


    def _load_primer_pairs(self, primer3_core_outfile):
        primer_pairs = {}

        f = pyfastaq.utils.open_file_read(primer3_core_outfile)
        while 1:
            results = self._get_next_primer3_sequence_results(f)
            if results is None:
                break

            if 'SEQUENCE_ID' not in results:
                raise Error('SEQUENCE_ID line not found in primer3_core output. Cannot continue')
            elif results['SEQUENCE_ID'] in primer_pairs:
                raise Error('Sequence name found twice:' + results['SEQUENCE_ID'] + ' ... cannot continue')
            elif self._primer3_found_no_primers(results):
                continue

            # remove everything after first white space so IDs match later,
            # as bowtie2 will do the same in its output SAM file
            results['SEQUENCE_ID'] = results['SEQUENCE_ID'].split(' ')[0]
            results_list = self._primer3_sequence_results_to_list(results)
            primer_pairs[results['SEQUENCE_ID']] = results_list

        pyfastaq.utils.close(f)

        for pair_list in primer_pairs.values():
            for pair in pair_list:
                if not pair.has_all_info():
                    raise Error('Error making primer pairs. Cannot continue')

        return primer_pairs


    def _write_primers_fasta(self, outfile):
        f = pyfastaq.utils.open_file_write(outfile)
        for genome_name in sorted(self.primer_pairs):
            for pair in self.primer_pairs[genome_name]:
                print(pair.left_fasta, file=f)
                print(pair.right_fasta, file=f)

        pyfastaq.utils.close(f)


    def run(self, outprefix):
        primer3_core_out = outprefix + '.primer3_core.out.gz'
        self._run_primer3_core(self.input_fasta, self.config_file, primer3_core_out)
        self.primer_pairs = self._load_primer_pairs(primer3_core_out)
        self._write_primers_fasta(outprefix + '.primers.fasta.gz')

