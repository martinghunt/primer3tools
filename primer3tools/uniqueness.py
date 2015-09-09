import os
import pyfastaq
import pysam
import primer3tools


class Error (Exception): pass


class PrimerUniqueness:
    def __init__(self, genomes_file, primer3_outdir, outprefix, min_product_length=50, max_product_length=1000):
        self.genomes_file = os.path.abspath(genomes_file)
        self.primer3_outdir = os.path.abspath(primer3_outdir)
        self.outprefix = outprefix
        self.min_product_length = min_product_length
        self.max_product_length = max_product_length


    def _cat_primer_fastas(self, genomes, outfile):
        original_line_length = pyfastaq.sequences.Fasta.line_length
        pyfastaq.sequences.Fasta.line_length = 0
        with open(outfile, 'w') as f_out:
           for genome_name in sorted(genomes):
               if genomes[genome_name].make_primers:
                   primers_fasta = os.path.join(self.primer3_outdir, genome_name + '.primers.fasta.gz')
                   file_reader = pyfastaq.sequences.file_reader(primers_fasta)
                   for seq in file_reader:
                       print(seq, file=f_out)

        pyfastaq.sequences.Fasta.line_length = original_line_length


    @staticmethod
    def _cat_all_genomes(genomes, outfile):
        with open(outfile, 'w') as f_out:
            for genome_name in sorted(genomes):
                file_reader = pyfastaq.sequences.file_reader(genomes[genome_name].fasta_file)
                for seq in file_reader:
                    seq.id = genome_name + '__' + seq.id
                    print(seq, file=f_out)


    @staticmethod
    def _is_perfect_hit(sam_record):
        if sam_record.is_unmapped:
            return False

        if sam_record.cigarstring != str(sam_record.query_length) + 'M':
            return False

        try:
            mismatches = sam_record.get_tag('NM')
        except:
            raise Error('No NM tag found for mapped read ' + sam_record.query_name)

        return mismatches == 0


    @staticmethod
    def _sam_to_read_sequence(sam_record):
        assert not sam_record.is_unmapped
        if sam_record.is_reverse:
            fa = pyfastaq.sequences.Fasta('x', sam_record.query_sequence)
            fa.revcomp()
            return fa.seq
        else:
            return sam_record.query_sequence


    @staticmethod
    def _filter_pairs_dict(pairs_dict):
        keys_to_remove = set()
        for name in pairs_dict:
            if 0 in [len(pairs_dict[name]['left']), len(pairs_dict[name]['right'])]:
                keys_to_remove.add(name)

        for name in keys_to_remove:
            del pairs_dict[name]


    @staticmethod
    def _parse_sam(infile):
        sam_reader = pysam.Samfile(infile, "r")
        pairs_dict = {}

        for read in sam_reader.fetch(until_eof=True):
            name_prefix, name_suffix = read.qname.rsplit('/', maxsplit=1)
            assert name_suffix in ['1', '2']
            left_or_right = 'left' if name_suffix == '1' else 'right'

            if name_prefix not in pairs_dict:
                pairs_dict[name_prefix] = {'left': {}, 'right': {}}

            if PrimerUniqueness._is_perfect_hit(read):
                refname = sam_reader.getrname(read.reference_id)
                if refname not in pairs_dict[name_prefix][left_or_right]:
                    pairs_dict[name_prefix][left_or_right][refname] = []
                pairs_dict[name_prefix][left_or_right][refname].append((read.reference_start, read.query_length, read.is_reverse, PrimerUniqueness._sam_to_read_sequence(read)))

        PrimerUniqueness._filter_pairs_dict(pairs_dict)
        return pairs_dict


    def _is_good_primer_pair(self, left_primer, right_primer):
        left_start, left_length, left_is_reverse, left_seq = left_primer
        right_start, right_length, right_is_reverse, right_seq = right_primer

        if left_is_reverse == right_is_reverse:
            return False

        if right_is_reverse:
            right_end = right_start + right_length - 1
            return left_start < right_start and self.min_product_length <= (right_end - left_start + 1) <= self.max_product_length
        else:
            left_end = left_start + left_length - 1
            return right_start < left_start and self.min_product_length <= (left_end - right_start + 1) <= self.max_product_length


    def _good_primer_pairs_from_lists(self, list1, list2):
        good_pairs = []

        for left_primer in list1:
            for right_primer in list2:
                if self._is_good_primer_pair(left_primer, right_primer):
                    good_pairs.append((left_primer, right_primer))

        return good_pairs


    def _all_primer_matches(self, hits_dict):
        matches = {}
        left_contig_names = set(hits_dict['left'].keys())
        right_contig_names = set(hits_dict['right'].keys())
        common_contig_names = left_contig_names.intersection(right_contig_names)

        for contig_name in common_contig_names:
            left_hits = hits_dict['left'][contig_name]
            right_hits = hits_dict['right'][contig_name]
            good_pairs = self._good_primer_pairs_from_lists(left_hits, right_hits)
            if len(good_pairs):
                matches[contig_name] = good_pairs

        return matches


    def _update_primer_hits(self, primer_hits, pairs_dict, genome_name):
        for primer_name_prefix in pairs_dict:
            matches = self._all_primer_matches(pairs_dict[primer_name_prefix])
            if len(matches):
                if primer_name_prefix not in primer_hits:
                    primer_hits[primer_name_prefix] = {}
                assert genome_name not in primer_hits[primer_name_prefix]
                primer_hits[primer_name_prefix][genome_name] = matches


    def _write_all_output_files(self, primer_hits, genomes):
        all_tsv = self.outprefix + '.all_primers.hits.tsv'
        unique_tsv = self.outprefix + '.unique_primers.tsv'
        f_out_all = pyfastaq.utils.open_file_write(all_tsv)
        f_out_unique = pyfastaq.utils.open_file_write(unique_tsv)
        genomes_with_unique_primer_pair = set()

        for primer_pair_prefix in sorted(primer_hits):
            all_hits = []
            unique_hit = None
            for genome_name in sorted(primer_hits[primer_pair_prefix]):
                for contig_name in sorted(primer_hits[primer_pair_prefix][genome_name]):
                    for left, right in primer_hits[primer_pair_prefix][genome_name][contig_name]:
                        unique_hit = [
                            genome_name,
                            contig_name,
                            str(left[0] + 1),
                            PrimerUniqueness._is_reverse_to_string(left[2]),
                            str(right[0] + 1),
                            PrimerUniqueness._is_reverse_to_string(right[2])
                        ]
                        all_hits.append(';'.join(unique_hit))
                        left_seq = left[3]
                        right_seq = right[3]

            if len(all_hits) == 1:
                print(unique_hit[0], unique_hit[1], unique_hit[2], unique_hit[4], left_seq, right_seq, primer_pair_prefix, sep='\t', file=f_out_unique)
                genomes_with_unique_primer_pair.add(genome_name)

            print(primer_pair_prefix, left_seq, right_seq, '\t'.join(all_hits), sep='\t', file=f_out_all)

        pyfastaq.utils.close(f_out_all)
        pyfastaq.utils.close(f_out_unique)

        f_out = pyfastaq.utils.open_file_write(self.outprefix + '.genome_uniqueness.tsv')
        for genome_name in sorted(genomes):
            if genomes[genome_name].make_primers:
                unique = '1' if genome_name in genomes_with_unique_primer_pair else '0'
                print(genome_name, unique, sep='\t', file=f_out)
        pyfastaq.utils.close(f_out)


    @staticmethod
    def _is_reverse_to_string(b):
        return '+' if not b else '-'


    def run(self):
        all_primers_fasta = self.outprefix + '.all_primers.fa'
        all_genomes_fasta = self.outprefix + '.all_genomes.fa'
        genomes = primer3tools.genome_set.GenomeSet(self.genomes_file)
        self._cat_primer_fastas(genomes, all_primers_fasta)
        primer_hits = {}

        for genome_name in genomes:
            genome = genomes[genome_name]
            bowtie_index = os.path.join(self.primer3_outdir, genome_name + '.bowtie2_index', 'index')
            sam_file = self.outprefix + '.tmp.' + genome_name + '.sam'
            primer3tools.mapping.run_bowtie2(all_primers_fasta, bowtie_index, sam_file)
            pairs_dict = self._parse_sam(sam_file)
            self._update_primer_hits(primer_hits, pairs_dict, genome_name)
            os.unlink(sam_file)

        self._write_all_output_files(primer_hits, genomes)
