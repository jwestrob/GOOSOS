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
    description='Given a directory of protein fasta files, extract the best hit for each HMM for each bin; get multifasta of each protein')
parser.add_argument('-fastadir', metavar='[PROTEIN FASTA DIRECTORY]',
                    help="Directory of protein fasta files to scan for ribosomal proteins. MUST END IN .faa EXTENSION")
parser.add_argument('-hmmdir', metavar='[HMM DIRECTORY]', help="Directory containing HMM files to scan with")
parser.add_argument('-outdir', metavar='[Output Directory]', default='output', help="Directory to store output files")
parser.add_argument('-evalue', metavar='[EVALUE THRESHOLD]', default=0.01,
        help="Evalue threshold for HMMsearch. Default: 0.01. (I pick the hit with the best e-value; \
        wait until alignment stage to curate)")
parser.add_argument('-threads', metavar='[NUM THREADS]', default=1, help="Number of threads to use")
parser.add_argument('-already_scanned', default=False, action='store_true', help='For if you already ran the HMMs')
parser.add_argument('-no_seqs', default=False, action='store_true', help='Dont pull out sequences to fasta')

def run_hmmsearch(protfile, hmmfile, wd, threshold):
    protein_id = protfile.split('/')[-1].split('.faa')[0]

    print('------------------------------------------------------------')
    print("Beginning HMMsearch...")
    print(protein_id, hmmfile)
    cmd = 'hmmsearch -o ' + wd + '/' + protein_id + '_' + hmmfile.split('/')[-1].split('.hmm')[0] + \
          '_hmmsearch.out  --notextw -E ' + str(threshold) + ' --cpu ' + str(1) + ' ' + hmmfile + \
          ' ' + protfile
    print(cmd)
    result = subprocess.getstatusoutput(cmd)
    if result[0] != 0:
        print('HMMsearch error (check for empty sequences in your protein FASTAs)')
        print('protein_id: ', protein_id)
        print('hmmfile: ', hmmfile)
        sys.exit()
    print(result)
    print('------------------------------------------------------------')
    return protein_id + '_' + hmmfile.split('/')[-1].split('.hmm')[0] + '_hmmsearch.out'


def extract_hits_by_outfile(dir, infile):
    hits = []
    e_values = []
    with open(dir + '/' + infile[0], 'r') as handle:
        for record in SearchIO.parse(handle, 'hmmer3-text'):
            hits.append(list(record))

    try:
        hits = hits[0]
    except:
        return

    good_hits = [hit._id for hit in hits]
    e_values = [hit.evalue for hit in hits]
    # If you have more than one hit, go with the hit that has the best e-value
    if len(good_hits) > 1:
        return good_hits[e_values.index(min(e_values))]
    else:
        try:
            return good_hits[0]
        except:
            return

def get_recs_for_fasta(hmm, fastadir):

    #Get name of FASTA so we can append that to the seqs for later identification
    fasta_id = fastadir.split('/')[-1]

    #There should only be one outfile matching the hmm provided
    hmmfile = list(filter(lambda x: hmm in x, os.listdir(fastadir)))[0]
    hits = []
    with open(fastadir + '/' + hmmfile, 'r') as handle:
        for record in SearchIO.parse(handle, 'hmmer3-text'):
            hits.append(list(record))

    try:
        hits = hits[0]
    except:
        return

    good_hits = [hit._id for hit in hits]

    out_recs = []

    fastafile = list(filter(lambda x: '.faa' in x, os.listdir(fastadir)))[0]
    fastafile = os.path.join(fastadir, fastafile)
    for rec in SeqIO.parse(fastafile, 'fasta'):
        if rec.id in good_hits:
            rec.id = fasta_id + '|' + rec.id
            out_recs.append(rec)

    return out_recs


def extract_hits_by_hmm(hmm, fastalist, outdir, threads):
    print("Extracting hits for " + hmm)
    p2 = Pool(threads)

    recs = list(p2.map(lambda fastaname: get_recs_for_fasta(hmm, outdir + '/' + fastaname), fastalist))

    return recs

def extract_hits(hmmlist, fastalist, outdir, threads):
    #List of recs (value to return)
    recs_by_hmm = []

    #I could use a map, but like... why
    for hmm in hmmlist:
        recs_by_hmm.append(extract_hits_by_hmm(hmm, fastalist, outdir, threads))

    return recs_by_hmm

def get_fastaheader_id(fasta):
    #Is this function necessary???
    for rec in SeqIO.parse(fasta, 'fasta'):
        if '.peg' in rec.id:
            id = rec.id.split('.peg')[0]
        #Everything with '.peg' will start with 'fig|', so this structure is necessary.
        #I usually just append the genome ID to the start of each header with a | delimiter before running.
        elif '|' in rec.id:
            id = rec.id.split('|')[0]
        else:
            print('Unrecognized header found. Aborting.')
            sys.exit()
        break
    return id

def make_hitstable_df(hits_by_hmm, hmmlist, fastalist, outdir):
    # Make matrix of zeros to store hits
    print("Making hits matrix...")
    hitstable = np.zeros((len(hmmlist), len(fastalist)))

    print("Len(hits_by_hmm):")
    print(len(hits_by_hmm))

    print("len(hmmlist): ", len(hmmlist))



    # Mark hits in table
    for hmm_idx, hmm in enumerate(hits_by_hmm):
        for genome_idx, genome_hits in enumerate(hmm):
            if type(genome_hits) is list:
                hits = len(genome_hits)
            #Used to make it a string if there was only one hit;
            #Not like that now but this doesn't hurt anything (all should be list)
            elif type(genome_hits) is str:
                hits = 1
            if genome_hits is None:
                hitstable[hmm_idx][genome_idx] = 0
            else:
                hitstable[hmm_idx][genome_idx] = hits


    hits = pd.DataFrame(hitstable).T
    hits.columns = hmmlist
    hits['id'] = fastalist

    #Get columns of DF
    cols = list(hits.columns.values)
    #Move IDs column to first index, for to make it look pretty
    cols.pop(cols.index('id'))
    hits = hits[['id'] + cols]
    #Write it to tsv in outdir without index (annoying)
    hits.to_csv(outdir + '/HITSTABLE.tsv', sep='\t', index=False)

def write_recs(recs_for_hmm, hmm_name, outdir):
    fasta_outdir = outdir + '/fastas'

    flatten = lambda l: [item for sublist in l for item in sublist]
    recs_for_hmm = flatten(recs_for_hmm)

    recs_for_hmm = list(filter(lambda x: type(x) is not None, recs_for_hmm))

    print("Writing recs for " + hmm_name)
    SeqIO.write(recs_for_hmm, fasta_outdir + '/' + hmm_name + '_hits.faa', 'fasta')
    return

def main():
    args = parser.parse_args()
    fastadir = str(Path(args.fastadir).absolute())
    hmmdir = str(Path(args.hmmdir).absolute())
    outdir = str(args.outdir)
    threshold = float(args.evalue)
    threads = int(args.threads)
    already_scanned = args.already_scanned
    no_seqs = args.no_seqs

    p = Pool(threads)

    # Make output directory
    if not os.path.exists(outdir):
        os.system('mkdir ' + outdir)
        outdir = str(Path(outdir).absolute())
    else:
        outdir = str(Path(outdir).absolute())

    # Get list of paths of all fastas
    fastalist_wpath = list(map(lambda file: os.path.join(fastadir, file), os.listdir(fastadir)))

    # Get list of all fastas
    fastalist = list(map(lambda file: file.split('.faa')[0], os.listdir(fastadir)))

    # Get list of paths of all HMM files
    hmmlist_wpath = list(map(lambda file: os.path.join(hmmdir, file), os.listdir(hmmdir)))

    # Get list of all HMMs
    hmmlist = list(map(lambda file: file.split('.hmm')[0], os.listdir(hmmdir)))

    hmm_outfiles = []

    # For each fasta, run all hmms
    if not already_scanned:
        for fastafile in fastalist_wpath:
            fastaoutdir = outdir + '/' + fastafile.split('/')[-1].split('.faa')[0]
            # Make outdir for HMMs
            if not os.path.exists(fastaoutdir):
                os.system('mkdir ' + outdir + '/' + fastafile.split('/')[-1].split('.faa')[0])
            #Make symbolic link
            os.system('ln -s ' + fastafile + ' ' + fastaoutdir + '/')
            hmm_outfiles.append([])

            # Run all HMMs for fastafile
            hmm_outfiles[-1] = list(p.map(lambda hmmfile: run_hmmsearch(fastafile, hmmfile, outdir, threshold), \
                                          hmmlist_wpath))

            # Move all outfiles to corresponding output directory
            for outfile in hmm_outfiles[-1]:
                os.system('mv ' + outdir + '/' + outfile + ' ' + fastaoutdir)

    # Make directory to store fastas
    if not os.path.exists(outdir + '/' + 'fastas'):
        os.system('mkdir ' + outdir + '/' + 'fastas')


    recs_list_by_hmm = extract_hits(hmmlist, fastalist, outdir, threads)

    make_hitstable_df(recs_list_by_hmm, hmmlist, fastalist, outdir)

    if not no_seqs:
        print("Getting recs and writing to fasta...")
        hmms_written = list(p.map(lambda hits:
                                        write_recs(hits,
                                        #Name of HMM to write (for fasta name)
                                        hmmlist[recs_list_by_hmm.index(hits)],
                                        outdir),
                                        #List of recs to iterate over
                                        recs_list_by_hmm))
        for hmm in hmmlist:
            if hmm not in hmms_written:
                print("Something went wrong. Hits weren't written for this HMM, please check:")
                print(hmm)

    print('boogie')


if __name__ == "__main__":
    main()
