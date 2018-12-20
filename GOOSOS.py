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
                    help="Directory of protein fasta files to scan for ribosomal proteins")
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


def get_recs_for_hits(hits_ids, hmm, fastadict, fastalist_wpath, fastalist, outdir):
    print("Hits IDs:")
    print(hits_ids)

    print('hmm:')
    print(hmm)

    sys.exit()
    hit_recs = []
    for hit in hits_ids:
        if hit is None:
            continue
        else:
            # print(hit)
            if '.peg' in hit:
                genome = fastadict[hit.split('.peg')[0]]
            else:
                genome = fastadict[hit.split('|')[0]]
            recs = list(SeqIO.parse(fastalist_wpath[fastalist.index(genome)], 'fasta'))
            hit_rec = list(filter(lambda x: x.id == hit, recs))[0]
            hit_rec.id = genome + '|' + hit
            hit_recs.append(hit_rec)
    print("Writing hits for: ", hmm)
    SeqIO.write(hit_recs, outdir + '/fastas/' + hmm + '.faa', 'fasta')
    return hmm

def get_hits_by_hmm(hmm, fastalist, outdir):

    #Make list of recs for this hmm (value to return)
    hmm_recs = []

    def get_hmm_by_fasta(fastaname, hmm):
        #Name of directory to find outfile/fasta in
        fastadir = outdir + '/' + fastaname

        #Get filename of relevant outfile
        hmmfile = list(filter(lambda x: hmm in x, os.listdir(fastadir)))[0]







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

    def get_fastaheader_id(fasta):
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

    #Get list of fasta header IDs by mapping to get_fastaheader_id fn
    fasta_header_ids = list(map(get_fastaheader_id, fastalist_wpath))

    #Make fasta dictionary (hopefully deprecated, let's see; dec 18 3:58 mountain time)
    fastadict = dict(zip(fasta_header_ids, fastalist))

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

    # Make matrix of zeros to store hits

    hits_by_hmm = []
    #Declare function to get hits for each HMM
    def extract_all_hits(fastaname, hmm):
        fastadir = outdir + '/' + fastaname
        #Get name of appropriate hmmfile, path
        hmmhits_for_fasta = list(filter(lambda x: hmm in x, os.listdir(fastadir)))
        hits = extract_hits_by_outfile(fastadir, hmmhits_for_fasta)
        return hits

    for hmm in hmmlist:
        print("Extracting hits for: ", hmm)
        relevant_outfiles = []
        hits_by_hmm.append([list(p.map(lambda fastaname:
                                        extract_all_hits(fastaname, hmm),
                                        fastalist)), hmm])


    print("Making hits matrix...")
    hitstable = np.zeros((len(hmmlist), len(fastalist)))

    # Mark hits in table
    for hmm_idx, hmm in enumerate(hits_by_hmm):
        for genome_idx, genome_hits in enumerate(hmm[0]):
            if type(genome_hits) is list:
                hits = len(genome_hits)
            elif type(genome_hits) is str:
                hits = 1
            if genome_hits is None:
                hitstable[hmm_idx][genome_idx] = 0
            else:
                hitstable[hmm_idx][genome_idx] = hits


    hits = pd.DataFrame(hitstable).T
    hits.columns = hmmlist
    hits['id'] = fastalist

    cols = list(hits.columns.values)
    cols.pop(cols.index('id'))
    hits = hits[['id'] + cols]
    hits.to_csv(outdir + '/HITSTABLE.tsv', sep='\t', index=False)

    if not no_seqs:
        print("Getting recs and writing to fasta...")
        hmms_written = list(p.map(lambda hits:
                   get_recs_for_hits(hits[0], hits[1], fastadict, fastalist_wpath, fastalist,
                                     outdir),
                   hits_by_hmm))
        for hmm in hmmlist:
            if hmm not in hmms_written:
                print(hmm)
        sys.exit()



    # recs_by_hmm = list(map(lambda hits: get_recs_for_hits(hits), hits_by_hmm))
    print('boogie')


if __name__ == "__main__":
    main()
