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
    aligns the sequences and concatenates them. Provides lots of histograms [outdir/figures]. \
    Also generates a partition file [outdir/partitions.nex] for use in phylogenomic analysis.')

parser.add_argument('-outdir', metavar='[GOOSOS output directory]', nargs=1,
                help="Provide the path to the directory you specified for GOOSOS output.")
parser.add_argument('-aln_name', metavar='[Name for concat alignment]', nargs=1,
                help="Provide a name for the concatenated alignment to generate. (Will be located in outdir/alignments).")
parser.add_argument('-exclude', metavar='[HMMs to exclude]', nargs='*', default=None,
                help="Any number of HMMs in your set that you'd like to exclude from the final alignment.")
parser.add_argument('-aln_concat', action='store_true', default=False,
                help="For if you already did filtering on your fastas (assumed to be in outdir/fastas) and just want to align/concatenate.")
parser.add_argument('-just_concat', action='store_true', default=False,
                help="For if you just want to concatenate some alignments (assumed to be in outdir/alignments). DOES NOT FILTER BASED ON THRESHOLD!")
parser.add_argument('-hits_threshold', metavar='[Lower threshold for num hits]', default=0.5,
                help="Percentage threshold (as a number between 0 and 1) \
                indicating how many hits out of the total a genome must have in \
                order to be included in the final alignment. Default: 50 (0.5)")
parser.add_argument('-inaccurate', action='store_true', default=False,
                help="Use faster (less accurate) parameters for MAFFT.")
parser.add_argument('-threads', metavar='[NUM THREADS]', default=1,
                help="Number of threads to use (helpful for MAFFT).")

def pass_sum_threshold(genome_id, threshold, df):
    total = len(df.columns.values)
    hits_threshold = float(threshold)*total
    num_hits = df[df.id == genome_id].sum(axis=1).tolist()[0]
    if num_hits < hits_threshold:
        return False
    else:
        return True

def align_fn(fastafile_wpath, outdir, threads, inaccurate):
    fastafile_id = fastafile_wpath.split('/')[-1].split('.faa')[0]
    if not inaccurate:
        #print('mafft --localpair --thread ' + str(threads) + ' --maxiterate 1000 ' + fastafile_wpath + ' > ' + outdir + '/alignments/' + fastafile_id + '_ALN.mfaa')
        os.system('mafft --localpair --reorder --thread ' + str(threads) + ' --maxiterate 1000 ' + fastafile_wpath + ' > '
                + outdir + '/alignments/' + fastafile_id + '_ALN.mfaa')
    else:
        os.system('mafft --reorder --thread ' + str(threads) + ' ' + fastafile_wpath + ' > '
                + outdir + '/' + alignments + '/' + fastafile_id + '_ALN.mfaa')
    return

def sort(alignment, genomes):
    alignment_length = len(alignment[0].seq)
    genomes_dict = dict(zip(genomes, range(len(genomes))))
    index_dict = dict(zip(range(len(genomes)), genomes))

    for rec in alignment:
        rec.id = rec.id.split('|')[0]

    out_alignment = []
    alignment_ids = [rec.id for rec in alignment]

    for genome in genomes:
        if genome in alignment_ids:
            #Wow this is convoluted but it's just finding the right index for this thing
            rec_to_append = alignment[alignment_ids.index(genome)]
            out_alignment.append(rec_to_append)
        else:
            #Alignment doesn't contain this genome; provide blank sequence
            blank_seq = '-'*alignment_length
            try:
                rec_to_append = SeqRecord(Seq(blank_seq), id=genome)
            except:
                print(genome)
            out_alignment.append(rec_to_append)
    return out_alignment

def concatenate(alignments_recs, alignment_name, alignments_dir):
    master_aln = alignments_recs[0]
    master_ids = [rec.id for rec in master_aln]
    for index, rec in enumerate(master_aln):
        for alignment in alignments_recs[1:]:
            master_aln[index].seq += alignment[index].seq

    for alignment in master_aln:
        alignment.description = ''

    fasta_extensions = ['.fa', '.faa', '.mfaa', '.fasta', '.fna']


    #Check to see if your filename ends with a .fa extension (or similar); if so just write the file as-is, else add extension
    if any(ext in alignment_name for ext in fasta_extensions):
        SeqIO.write(master_aln, alignments_dir + '/' + alignment_name, 'fasta')
    else:
        SeqIO.write(master_aln, alignments_dir + '/' + alignment_name + '.mfaa', 'fasta')
    return

def make_partition_file(alignments_recs, alignments_dir):
    outfile_list = []
    outfile_list.append('#nexus\n')
    #print('#nexus')
    outfile_list.append('begin sets;\n')
    #print('begin sets;')
    position = 1
    for aln_num, alignment in enumerate(alignments_recs):
        line = '\tcharset part' + str(aln_num) +  ' = ' + str(int(position)) + '-' + str(int(position) + len(alignment[0].seq)) + ';\n'
        #print(line)
        position = int(position) + len(alignment[0].seq)
        outfile_list.append(line)
    outfile_list.append('end;')
    with open(alignments_dir + '/partitions.nex', 'w') as outfile:
        outfile.writelines(outfile_list)


def throw_flags(hitstable, genomes_passed_threshold):
    num_hits_by_hmm = hitstable.sum(axis=0)
    threshold = float(0.25)*float(len(hitstable))
    for element in num_hits_by_hmm:
        if float(element) < threshold:
            print('abawaca4')
    #DO THIS LATER
    return





def main(args):
    outdir = str(Path(args.outdir[0]).absolute())
    aln_concat = args.aln_concat
    just_concat = args.just_concat
    threshold = args.hits_threshold
    exclude = args.exclude
    threads = args.threads
    inaccurate = args.inaccurate
    alignment_name = args.aln_name[0]

    if threshold is not None and just_concat == True:
        print("Since you specified the -just_concat flag, the threshold will not be used to filter.")
        print("This is assumed to have been done already in a previous run.")

    if aln_concat and just_concat:
        print("You can't use the -aln_concat and -justconcat flags at the same time. See the README.")
        sys.exit()

    if not aln_concat and not just_concat:

        fastas = os.listdir(outdir + '/fastas')
        #If you have aln_concat flag activated, this didn't happen earlier
        fastas_wpath = list(map(lambda x: os.path.join(outdir + '/fastas', x), fastas))

        if exclude is not None:
            fastas_wpath = list(filter(lambda x: x.split('/')[-1].split('.faa')[0] not in exclude,
                                                 fastas_wpath))

        #Get df of all hits, just in case
        all_df = pd.read_csv(outdir + '/all_hits_evalues_df.tsv', sep='\t')

        #Get hitstable
        hitstable = pd.read_csv(outdir + '/HITSTABLE.tsv', sep='\t')
        hitstable.columns.values[0] = 'id'

        #Get order of genomes
        genomes = hitstable['id'].tolist()


        genomes_passed_threshold = list(filter(lambda genome_id: pass_sum_threshold(genome_id, threshold, hitstable),
                                                              genomes))
        if outdir.endswith('/'):
            fastadir = outdir + 'fastas'
        else:
            fastadir = outdir + '/fastas'

        #You don't have to be loonelyyyyy at fastasonly.com
        fastas_only = list(filter(lambda x: '.faa' in x, os.listdir(fastadir)))
        fastas_recs = list(map(lambda fastafile:
                      list(SeqIO.parse(os.path.join(fastadir, fastafile),'fasta')),
                       fastas_only))


        orf_id_list = all_df.orf_id.tolist()

        fastas_recs_filtered = []
        for fasta in fastas_recs:
            fasta = list(filter(lambda x: x.id.split('|')[0] in genomes_passed_threshold, fasta))
            newrecs = []
            for rec in fasta:
                if rec.id.split('|')[-1] in orf_id_list:
                    newrecs.append(rec)
            fastas_recs_filtered.append(newrecs)

        if sum([len(x) for x in fastas_recs_filtered]) != len(orf_id_list):
            print("WHOOPS")
            sys.exit()
        else:
            print("OK!")

        #Make subdirectory to store filtered fastas
        goodseqs = fastadir + '/filtered_fastas'
        if not os.path.exists(goodseqs):
            os.mkdir(goodseqs)

        for index, fastafile in enumerate(fastas_only):
            SeqIO.write(fastas_recs_filtered[index], os.path.join(goodseqs,fastafile), 'fasta')

        #throw_flags(hitstable, genomes_passed_threshold)

        print(str(len(genomes_passed_threshold)) + " genomes passed the threshold for number of hits.")


    if not just_concat:
        if aln_concat:
            fastas = os.listdir(outdir + '/fastas')
            #If you have aln_concat flag activated, this didn't happen earlier
            fastas_wpath = list(map(lambda x: os.path.join(outdir + '/fastas', x), fastas))

            if exclude is not None:
                fastas_wpath = list(filter(lambda x: x.split('/')[-1].split('.faa')[0] not in exclude,
                                                    fastas_wpath))
            print("Beginning alignments...")
            if not os.path.exists(outdir + '/alignments'):
                os.system('mkdir ' + outdir + '/alignments')
            list(map(lambda x: align_fn(x, outdir, threads, inaccurate), fastas_wpath))
        else:
            filtered_fastas = os.listdir(goodseqs)
            filtered_fastas_wpath = list(map(lambda x: os.path.join(goodseqs, x), filtered_fastas))
            print("Beginning alignments...")
            if not os.path.exists(outdir + '/alignments'):
                os.system('mkdir ' + outdir + '/alignments')
            list(map(lambda x: align_fn(x, outdir, threads, inaccurate), filtered_fastas_wpath))

    alignments_dir = outdir + '/alignments'
    alignments = list(filter(lambda x: 'ALN' in x, os.listdir(alignments_dir)))
    alignments_recs = list(map(lambda alignment:
                        list(SeqIO.parse(os.path.join(alignments_dir, alignment), 'fasta')),
                        alignments))



    #Does this still work properly if you skip steps?? Check this later
    #Of course it doesn't lmao good job past jacob

    if not just_concat:
        alignments_recs_sorted = list(map(lambda alignment: sort(alignment, genomes_passed_threshold),
                                                            alignments_recs))
    else:
        #It's assumed that you've filtered at this point
        redundant_ids = []
        for alignment in alignments_recs:
            for rec in alignment:
                redundant_ids.append(rec.id.split('|')[0])

        all_ids = list(set(redundant_ids))

        alignments_recs_sorted = list(map(lambda alignment: sort(alignment, all_ids),
                                                            alignments_recs))

    lengths = [len(alignment[0].seq) for index, alignment in enumerate(alignments_recs_sorted)]
    print("Your concatenated alignment is " + str(sum(lengths)) + ' characters long. Congratulations')

    #Make partition file for phylogenomic analysis
    make_partition_file(alignments_recs_sorted, alignments_dir)

    #Concatenate and make final alignment
    concatenate(alignments_recs_sorted, alignment_name, alignments_dir)



if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
