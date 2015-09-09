# primer3tools

This is a tool to design primers for a set of genomes, and check the uniqueness
of the primers against a second (presumably larger) set of background genomes.
Primer sequences are generated using Primer3, and uniqueness is checked by
mapping with Bowtie2.


# Installation

primer3tools has the following dependencies:
  * [primer3] [Primer3] version >= 2.3.6 (`primer3_core` must be in your path)
  * [bowtie2] [Bowtie2] version >= 2.1.0


Download the latest release from this github repository,
or clone the repository. Then run the tests:

    python3 setup.py test

If the tests all pass, install:

    python3 setup.py install


Or install to a location of your choice with:

    python3 setup.py install --prefix /path/to/installation/dir/



# Usage

The installation installs a script called `primer3tools`. The workflow is as follows.


## Make a primer3 config file

Make a config file called `primer3.config` by running:

    primer3tools make_config primer3.config

The fourth line of this file, which looks like this:

    PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/absolute/path/to/src/primer3_config/

needs editing to point to a valid path. To use primer3 default settings, replace
`/absolute/path/to/` with your local primer3 installation directory.

For internal Sanger users only, use this:

    PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/software/pathogen/external/apps/usr/local/primer3-2.3.6/src/primer3_config/


## Run primer3 (and bowtie2 indexing) on a set of genomes

Make a config file that has information on the genomes. The format is one genome per line, with three tab-delimited columns:

1. Name of the genome, which must be unique
2. Name of FASTA file
3. A 0 or a 1. 1 means that primers are to be designed for this genome.
   0 means that primers are not designed, it is just used when primers for
   other genomes are checked for uniqueness later.

Run with:

    primer3tools batch primer3.config genomes.config Output_directory

Multiple threads can be used using the option `--threads` (default is 1), for example to use 4 threads:

    primer3tools batch --threads 4  primer3.config genomes.config Batch_output_directory


## Check uniqueness of primers

This uses the output directory made by the previous stage to check the uniqueness of all the primers
reported by primer3. Run it with:

    primer3tools get_unique genomes.config Batch_output_directory out

The minimum and maximum allowed PCR product length can be changed using the options
`--min_product_length` and `--max_product_length`.

The output files are called `out.*`. These are:

* **`out.all_primers.fa`** - a FASTA file of all the primer pairs reported by primer3. The name of each
  primer is: `genome-name__contig-name__primer3-number__fwd-start__rev-start/1_or_2`, where:
  * `genome-name` is the name of the genome
  * `contig-name` is the name of the contig
  * `primer3-number` is the primer pair number for this contig, given by primer3
  * `fwd-start` is the start position of the forwards primer
  * `rev-start` is the start position of the reverse primer
  * `1_or_2` - `1` means the forward primer and `2` means the reverse primer

* **`out.all_primers.hits.tsv`** - a tab-separated file showing all the matches each primer had to
  any of the input genomes. Currently, a match is defined by bowtie2 reporting a perfect match, and
  the forward and reverse sequences distance and orientation must be correct. There is
  one line per primer pair. The columns are:
  1. Primer pair name (see description of name in `out.all_primers.fa`).
  2. Forward primer sequence
  3. Reverse primer sequence
  4. Position of match to genome (as reported by bowtie2).
     This is in the form: genome-name;contig-name;fwd-start;forward-strand;rev-start;reverse-strand.

  Primer pairs with more than one hit will have an extra column per hit, in the same format as column 4.

* **`out.unique_primers.tsv`** - a tab-delimited file of primer pairs that have exactly one match
across all input genomes. The columns are:
  1. Genome name
  2. Contig name
  3. Start position of forwards primer
  4. Start position of reverse primer
  5. Forward primer sequence
  6. Reverse primer sequence
  7. Primer pair name (matches name used in `out.all_primers.fa` and `out.all_primers.hits.tsv`)

* **`out.genome_uniqueness.tsv`** - a tab-delimited file showing which genomes had at least one
  unique primer pair. It lists all of the genome names that had primer3 run on them (column 1). The
  second column has either a 1 or a 0. 1 means at least one unique primer pair was found, otherwise
  the second column has 0.


  [bowtie2]: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
  [primer3]: http://sourceforge.net/projects/primer3/
