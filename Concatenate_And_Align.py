from pathos.multiprocessing import ProcessingPool as Pool
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, SearchIO
from Bio.Alphabet import IUPAC
import os, sys, pandas as pd
from pathlib import Path
from Bio.Seq import Seq
import numpy as np
import subprocess
import argparse

parser = argparse.ArgumentParser(
    description='So you ran GOOSOS and got some sorted fasta files. \
    This script filters out genomes and genes that are below inclusion thresholds, \
    aligns the sequences and concatenates them. Provides lots of histograms (outdir/figures).')

parser.add_argument('-outdir', metavar='[GOOSOS output directory]', nargs=1,
                help="Provide the path to the directory you specified for GOOSOS output.")
parser.add_argument('-exclude', metavar='[HMMs to exclude]', nargs='*',
                help="Any number of HMMs in your set that you'd like to exclude from the final alignment.")
parser.add_argument('-aln_concat', action='store_true', default=False,
                help="For if you already did filtering on your fastas (assumed to be in outdir/fastas) and just want to align/concatenate.")
parser.add_argument('-justconcat', action='store_true', default=False,
                help="For if you just want to concatenate some alignments (assumed to be in outdir/alignments).")
parser.add_argument('-hits_threshold', metavar='[Lower threshold for num hits]', default=0.5,
                help="Percentage threshold (as a number between 0 and 1) \
                indicating how many hits out of the total a genome must have in \
                order to be included in the final alignment. Default: 50% (0.5)")


args = parser.parse_args()

outdir = args.outdir
aln_concat = args.aln_concat
just_concat = args.justconcat
hits_threshold = args.hits_threshold

if aln_concat and just_concat:
    print("You can't use the -aln_concat and -justconcat flags at the same time. See the README.")
    sys.exit()

all_df = pd.read_csv(outdir + '/all_hits_evalues_df.tsv', sep='\t')
hitstable = pd.read_csv(outdir + '/HITSTABLE.tsv', sep='\t')
