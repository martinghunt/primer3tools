import argparse
import primer3tools

def run():
    parser = argparse.ArgumentParser(
        description = 'Make default primer3_core config file',
        usage = 'primer3tools <outfile>',
        epilog = 'IMPORTANT: output file must be edited to fix PRIMER_THERMODYNAMIC_PARAMETERS_PATH option',
    )

    parser.add_argument('outfile', help='Name of output file')
    options = parser.parse_args()

    with open(options.outfile, 'w') as f:
        print(r'''Primer3 File - http://primer3.sourceforge.net
P3_FILE_TYPE=settings

PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/absolute/path/to/src/primer3_config/
PRIMER_FIRST_BASE_INDEX=1
PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1
PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=0
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_LIBERAL_BASE=1
PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0
PRIMER_LOWERCASE_MASKING=0
PRIMER_PICK_ANYWAY=1
PRIMER_EXPLAIN_FLAG=1
PRIMER_TASK=generic
PRIMER_MIN_QUALITY=0
PRIMER_MIN_END_QUALITY=0
PRIMER_QUALITY_RANGE_MIN=0
PRIMER_QUALITY_RANGE_MAX=100
PRIMER_MIN_SIZE=18
PRIMER_OPT_SIZE=20
PRIMER_MAX_SIZE=23
PRIMER_MIN_TM=57.0
PRIMER_OPT_TM=59.0
PRIMER_MAX_TM=62.0
PRIMER_PAIR_MAX_DIFF_TM=5.0
PRIMER_TM_FORMULA=1
PRIMER_PRODUCT_MIN_TM=-1000000.0
PRIMER_PRODUCT_OPT_TM=0.0
PRIMER_PRODUCT_MAX_TM=1000000.0
PRIMER_MIN_GC=30.0
PRIMER_OPT_GC_PERCENT=50.0
PRIMER_MAX_GC=70.0
PRIMER_PRODUCT_SIZE_RANGE=150-250 100-300 301-400 401-500 501-600 601-700 701-850 851-1000
PRIMER_NUM_RETURN=5
PRIMER_MAX_END_STABILITY=9.0
PRIMER_MAX_LIBRARY_MISPRIMING=12.00
PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=20.00
PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00
PRIMER_MAX_SELF_ANY_TH=45.0
PRIMER_MAX_SELF_END_TH=35.0
PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0
PRIMER_PAIR_MAX_COMPL_END_TH=35.0
PRIMER_MAX_HAIRPIN_TH=24.0
PRIMER_MAX_TEMPLATE_MISPRIMING=12.00
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00
PRIMER_MAX_SELF_ANY=8.00
PRIMER_MAX_SELF_END=3.00
PRIMER_PAIR_MAX_COMPL_ANY=8.00
PRIMER_PAIR_MAX_COMPL_END=3.00
PRIMER_MAX_NS_ACCEPTED=0
PRIMER_MAX_POLY_X=4
PRIMER_INSIDE_PENALTY=-1.0
PRIMER_OUTSIDE_PENALTY=0
PRIMER_GC_CLAMP=0
PRIMER_MAX_END_GC=5
PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE=3
PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE=3
PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION=7
PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION=4
PRIMER_SALT_MONOVALENT=50.0
PRIMER_SALT_CORRECTIONS=1
PRIMER_SALT_DIVALENT=1.5
PRIMER_DNTP_CONC=0.6
PRIMER_DNA_CONC=50.0
PRIMER_SEQUENCING_SPACING=500
PRIMER_SEQUENCING_INTERVAL=250
PRIMER_SEQUENCING_LEAD=50
PRIMER_SEQUENCING_ACCURACY=20
PRIMER_WT_SIZE_LT=1.0
PRIMER_WT_SIZE_GT=1.0
PRIMER_WT_TM_LT=1.0
PRIMER_WT_TM_GT=1.0
PRIMER_WT_GC_PERCENT_LT=0.0
PRIMER_WT_GC_PERCENT_GT=0.0
PRIMER_WT_SELF_ANY_TH=0.0
PRIMER_WT_SELF_END_TH=0.0
PRIMER_WT_HAIRPIN_TH=0.0
PRIMER_WT_TEMPLATE_MISPRIMING_TH=0.0
PRIMER_WT_SELF_ANY=0.0
PRIMER_WT_SELF_END=0.0
PRIMER_WT_TEMPLATE_MISPRIMING=0.0
PRIMER_WT_NUM_NS=0.0
PRIMER_WT_LIBRARY_MISPRIMING=0.0
PRIMER_WT_SEQ_QUAL=0.0
PRIMER_WT_END_QUAL=0.0
PRIMER_WT_POS_PENALTY=0.0
PRIMER_WT_END_STABILITY=0.0
PRIMER_PAIR_WT_PRODUCT_SIZE_LT=0.0
PRIMER_PAIR_WT_PRODUCT_SIZE_GT=0.0
PRIMER_PAIR_WT_PRODUCT_TM_LT=0.0
PRIMER_PAIR_WT_PRODUCT_TM_GT=0.0
PRIMER_PAIR_WT_COMPL_ANY_TH=0.0
PRIMER_PAIR_WT_COMPL_END_TH=0.0
PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH=0.0
PRIMER_PAIR_WT_COMPL_ANY=0.0
PRIMER_PAIR_WT_COMPL_END=0.0
PRIMER_PAIR_WT_TEMPLATE_MISPRIMING=0.0
PRIMER_PAIR_WT_DIFF_TM=0.0
PRIMER_PAIR_WT_LIBRARY_MISPRIMING=0.0
PRIMER_PAIR_WT_PR_PENALTY=1.0
PRIMER_PAIR_WT_IO_PENALTY=0.0
PRIMER_INTERNAL_MIN_SIZE=18
PRIMER_INTERNAL_OPT_SIZE=20
PRIMER_INTERNAL_MAX_SIZE=27
PRIMER_INTERNAL_MIN_TM=57.0
PRIMER_INTERNAL_OPT_TM=60.0
PRIMER_INTERNAL_MAX_TM=63.0
PRIMER_INTERNAL_MIN_GC=20.0
PRIMER_INTERNAL_OPT_GC_PERCENT=50.0
PRIMER_INTERNAL_MAX_GC=80.0
PRIMER_INTERNAL_MAX_SELF_ANY_TH=47.00
PRIMER_INTERNAL_MAX_SELF_END_TH=47.00
PRIMER_INTERNAL_MAX_HAIRPIN_TH=47.00
PRIMER_INTERNAL_MAX_SELF_ANY=12.00
PRIMER_INTERNAL_MAX_SELF_END=12.00
PRIMER_INTERNAL_MIN_QUALITY=0
PRIMER_INTERNAL_MAX_NS_ACCEPTED=0
PRIMER_INTERNAL_MAX_POLY_X=5
PRIMER_INTERNAL_MAX_LIBRARY_MISHYB=12.00
PRIMER_INTERNAL_SALT_MONOVALENT=50.0
PRIMER_INTERNAL_DNA_CONC=50.0
PRIMER_INTERNAL_SALT_DIVALENT=1.5
PRIMER_INTERNAL_DNTP_CONC=0.0
PRIMER_INTERNAL_WT_SIZE_LT=1.0
PRIMER_INTERNAL_WT_SIZE_GT=1.0
PRIMER_INTERNAL_WT_TM_LT=1.0
PRIMER_INTERNAL_WT_TM_GT=1.0
PRIMER_INTERNAL_WT_GC_PERCENT_LT=0.0
PRIMER_INTERNAL_WT_GC_PERCENT_GT=0.0
PRIMER_INTERNAL_WT_SELF_ANY_TH=0.0
PRIMER_INTERNAL_WT_SELF_END_TH=0.0
PRIMER_INTERNAL_WT_HAIRPIN_TH=0.0
PRIMER_INTERNAL_WT_SELF_ANY=0.0
PRIMER_INTERNAL_WT_SELF_END=0.0
PRIMER_INTERNAL_WT_NUM_NS=0.0
PRIMER_INTERNAL_WT_LIBRARY_MISHYB=0.0
PRIMER_INTERNAL_WT_SEQ_QUAL=0.0
PRIMER_INTERNAL_WT_END_QUAL=0.0
=''', file=f)

