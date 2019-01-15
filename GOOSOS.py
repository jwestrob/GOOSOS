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
    description='Given a directory of protein fasta files, extract the hits for each HMM for each bin; get multifasta of each protein')

parser.add_argument('-nuc', action='store_true',
                    default=False, help='Run GOOSOS Nucleotide workflow. Ignores -protdir')
parser.add_argument('-prot', action='store_true',
                    default=False, help='Run GOOSOS protein  workflow. If -nucdir is provided, includes predicted proteins from nucleotide sequences.')
parser.add_argument('-nucdir', metavar='[NUCLEOTIDE FASTA DIRECTORY]',
                    help='Directory containing nucleotide sequences.')
parser.add_argument('-protdir', metavar='[PROTEIN FASTA DIRECTORY]',
                    help="Directory of protein fasta files to scan for ribosomal proteins. MUST END IN .faa EXTENSION")
parser.add_argument('-hmmdir', metavar='[HMM DIRECTORY]', help="Directory containing HMM files to scan with")
parser.add_argument('-outdir', metavar='[Output Directory]', default='output', help="Directory to store output files")
parser.add_argument('-evalue', metavar='[EVALUE THRESHOLD]', default=0.01,
                    help="Evalue threshold for HMMsearch. Default: 0.01.")
parser.add_argument('-threads', metavar='[NUM THREADS]', default=1, help="Number of threads to use")
parser.add_argument('-already_predicted', default=False, action='store_true', help='For if you already ran Prodigal.')
parser.add_argument('-already_scanned', default=False, action='store_true', help='For if you already ran the HMMs')
parser.add_argument('-no_seqs', default=False, action='store_true', help='Dont pull out sequences to fasta')
parser.add_argument('-best', default=False, action='store_true', help='Only pull out best hit per genome')

def hmmpress(hmmlist_wpath, outdir):
    #Concatenate all hmm files together and press them into a binary
    list_of_hmms = ' '.join(hmmlist_wpath)

    #Make folder to store hmmpress files in
    os.mkdir(outdir + '/hmmpress')
    cwd = os.getcwd()

    os.chdir(outdir)
    os.system('cat ' + list_of_hmms + ' > concatenated_hmms.hmm')

    os.system('hmmpress concatenated_hmms.hmm')
    os.system('mv ' + outdir + '/concatenated_hmms.* ' + outdir + '/hmmpress/')

    os.chdir(cwd)
    return

def run_hmmsearch(protfile, hmmfile, outdir, threshold):
    protein_id = protfile.split('/')[-1].split('.faa')[0]

    #print(protein_id, hmmfile)
    cmd = 'hmmsearch -o ' + outdir + '/hmmsearch/' + protein_id + '_' + hmmfile.split('/')[-1].split('.hmm')[0] + \
          '_hmmsearch.out  --notextw -E ' + str(threshold) + ' --cpu ' + str(1) + ' ' + hmmfile + \
          ' ' + protfile
    #print(cmd)
    result = subprocess.getstatusoutput(cmd)
    if result[0] != 0:
        print('HMMsearch error (check for empty sequences in your protein FASTAs)')
        print('protein_id: ', protein_id)
        #print('hmmfile: ', hmmfile)
        sys.exit()
    #print(result)
    return protein_id + '_' + hmmfile.split('/')[-1].split('.hmm')[0] + '_hmmsearch.out'

def run_hmmscan(protfile, outdir, threshold):
    genome_id = protfile.split('/')[-1].split('.faa')[0].split('.fna')[0].split('.fa')[0].split('.fasta')[0]

    #print(protein_id, hmmfile)
    cmd = 'hmmscan --domtblout ' + outdir + '/hmmscan/' + genome_id + '/' + genome_id + '_hmmsearch.out  --notextw -E ' \
            + str(threshold) + ' --cpu ' + str(1) + ' ' + outdir + '/hmmpress/concatenated_hmms.hmm ' + protfile + ' > /dev/null 2>&1'
    #print(cmd)
    result = subprocess.run(cmd, shell=True, check=True)
    if result.returncode != 0:
        print(result)
        print('HMMscan error (check for empty sequences in your protein FASTAs)')
        print('genome_id: ', genome_id)
        #print('hmmfile: ', hmmfile)
        sys.exit()
    #print(result)
    #Parse file with awk/perl nonsense; generate .parse file
    return parse_hmmdomtbl(outdir, genome_id + '_hmmsearch.out', threshold)


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


def get_rec_for_hit(genome_id, orf, outdir):
    genome_dir = outdir + '/hmmscan/' + genome_id + '/'
    protfile = list(filter(lambda x: '.faa' in x, os.listdir(genome_dir)))[0]
    genome_recs = list(SeqIO.parse(genome_dir + protfile, 'fasta'))

    desired_hit = list(filter(lambda x: orf in x.id, genome_recs))[0]
    desired_hit.id = genome_id + '|' + desired_hit.id

    return desired_hit

def extract_hits_by_hmm(red_df, threads, outdir):
    print("Extracting hits for " + red_df.iloc[0].family_hmm)
    p2 = Pool(threads)

    #Make list of genome_id / orf pairs
    id_orf_list = red_df[['genome_id', 'orf_id']].values.tolist()


    recs = list(p2.map(lambda id_orf: get_rec_for_hit(id_orf[0], id_orf[1], outdir), id_orf_list))

    return recs

def extract_hits(all_df, threads, outdir):
    #List of recs (value to return)
    recs_by_hmm = []

    #I could use a map, but like... why
    for hmm in all_df['family_hmm'].unique().tolist():
        red_df = all_df[all_df['family_hmm'] == hmm]
        recs_by_hmm.append([extract_hits_by_hmm(red_df, threads, outdir), hmm])

    return recs_by_hmm

def make_hitstable_df(recs_by_hmm, hmmlist, fastalist, outdir):

    # Make matrix of zeros to store hits
    print("Making hits matrix...")
    hitstable = np.zeros((len(hmmlist), len(fastalist)))

    hits = pd.DataFrame(hitstable).T
    hits.columns = hmmlist
    print(fastalist)
    sys.exit()
    hits['id'] = fastalist
    #Get columns of DF
    cols = list(hits.columns.values)
    #Move IDs column to first index, for to make it look pretty
    cols.pop(cols.index('id'))
    hits = hits[['id'] + cols]

    print(hits)
    sys.exit()

    # Mark hits in table
    for hmm_recs, hmm in recs_by_hmm:

        hmm_idx = hmmlist.index(hmm)

        for genome_hit in hmm_recs:
            #Extract genome ID from fasta header
            try:
                genome_id = genome_hit.id.split('|')[0]
            except:
                print(genome_id)
                print(type(genome_id))
            print(hits[hits['id'] == genome_id])
            hits[hits['id'] == genome_id][hmm] += 1


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

def run_prodigal(fastafile_wpath, outdir):
    #print('prodigal -i '+ fastafile_wpath + ' -a ' + outdir + '/proteins/' +
    #                        fastafile_wpath.split('/')[-1] + '.faa -m -p single > /dev/null 2>&1')
    os.system('prodigal -i '+ fastafile_wpath + ' -a ' + outdir + '/proteins/' +
                            fastafile_wpath.split('/')[-1] + '.faa -m -p single > /dev/null 2>&1')
    print('Genes predicted for ' + fastafile_wpath.split('/')[-1])
    return

def parse_hmmdomtbl(outdir, hmmoutfile, threshold):

    """
    Takes output file from HMMscan (--domtblout), parses it, and yields a
    genome-specific dataframe of each hit (above the evalue threshold) for each
    domain.
    """

    genome_id = hmmoutfile.split('_hmmsearch.out')[0].split('.fasta')[0].split('.fna')[0].split('.fa')[0]
    hmmoutfile_wpath = outdir + '/hmmscan/' + genome_id + '/' + hmmoutfile

    with open(hmmoutfile_wpath, 'r') as infile:
        lines = infile.readlines()

    domtbl_header = ['target name', 'accession', 'tlen', 'query name', 'accession', 'qlen', 'E-value',
                    'full_score',  'full_bias',   'dom_#',  'dom_of',  'c-Evalue',  'i-Evalue',  'dom_score',  'dom_bias',  'hmm_from',
                    'hmm_to',  'ali_from', 'ali_to',  'env_from', 'env_to',  'mean_posterior', 'description of target']

    desired_header = ['family_hmm', 'hmm_length', 'orf_id',
                      'query_length', 'evalue', 'hmm_start',
                      'hmm_end', 'query_start', 'query_end']

    #Remove all lines with '#' beginning character
    lines_filtered = list(filter(lambda x: x[0] != '#', lines))
    if len(lines_filtered) == 0:
        print("No hits for " + genome_id)
        orflist_header = ['family_hmm', 'hmm_length', 'orf_id', 'overall_evalue', 'dom1_cevalue', 'dom1_hmmstart',
                          'dom1_hmmend', 'dom1_querystart', 'dom1_queryend', 'dom2_evalue', 'dom2_hmmstart',
                          'dom2_hmmend', 'dom2_querystart', 'dom2_queryend']
        empty_df = pd.DataFrame(columns = orflist_header)
        empty_df.to_csv(outdir + '/hmmscan/' + genome_id + '/' + genome_id + '.parse', sep='\t', index=False)
        return outdir + '/hmmscan/' + genome_id + '/' + genome_id + '.parse'
    #Remove newline characters and split by whitespace
    lines_filtered = list(map(lambda x: x.strip('\n').split(), lines_filtered))
    #Python, by default, splits the description at the end, causing dimension mismatch.
    #Let's get rid of the extra elements.
    lines_filtered = list(map(lambda x: x[0:23], lines_filtered))
    #print(domtbl_header)
    #print(lines_filtered[0])
    #sys.exit()
    #Make pandas DF to store lines, then add column names
    lines_df = pd.DataFrame(lines_filtered, columns=domtbl_header)


    #Make DF to store properly arranged data
    goodheader_df = pd.DataFrame(columns=desired_header)

    goodheader_df['family_hmm'] = lines_df['target name']
    #Insert genome ID to avoid confusion
    goodheader_df['genome_id'] = pd.Series([genome_id]*len(lines_df))
    goodheader_df['hmm_length'] = lines_df['tlen']
    goodheader_df['orf_id'] = lines_df['query name']
    goodheader_df['query_length'] = lines_df['qlen']
    goodheader_df['evalue'] = lines_df['E-value']
    goodheader_df['c_evalue'] = lines_df['c-Evalue']
    goodheader_df['hmm_start'] = lines_df['hmm_from']
    goodheader_df['hmm_end'] = lines_df['hmm_to']
    goodheader_df['query_start'] = lines_df['ali_from']
    goodheader_df['query_end'] = lines_df['ali_to']

    unique_orfs = goodheader_df['orf_id'].unique()
    #print("Unique orfs: ")
    #print(unique_orfs)

    orflist = []
    orflist_header = ['family_hmm', 'genome_id', 'orf_id', 'hmm_length', 'overall_evalue', 'dom1_cevalue', 'dom1_hmmstart',
                      'dom1_hmmend', 'dom1_querystart', 'dom1_queryend', 'dom2_evalue', 'dom2_hmmstart',
                      'dom2_hmmend', 'dom2_querystart', 'dom2_queryend']
    for orf in unique_orfs:
        red_df = goodheader_df[goodheader_df['orf_id'] == orf]

        #For the unlikely scenario where you have more than one HMM hit on a single ORF
        if len(red_df['family_hmm'].unique()) > 1:
            #You want ascending=True because you want the smallest values first
            sorted_red = red_df.sort_values(by='evalue', ascending=True)
            best_ORF = sorted_red.iloc[0].family_hmm
            red_df = red_df[red_df['family_hmm'] == best_ORF]

        # if len(red_df) > 2: Complain to me about it and I can do something later
        # or, I REALLY DON'T CARE DO U?
        if len(red_df) >= 2:
            sorted_red = red_df.sort_values(by='evalue', ascending=True)
            goodrow = sorted_red.iloc[0]
            worse_row = sorted_red.iloc[1]
            if float(goodrow.evalue) < threshold:
                orflist.append([goodrow.family_hmm,
                                genome_id,
                                goodrow.orf_id,
                                goodrow.hmm_length,
                                goodrow.evalue,
                                goodrow.c_evalue,
                                goodrow.hmm_start,
                                goodrow.hmm_end,
                                goodrow.query_start,
                                goodrow.query_end,
                                worse_row.c_evalue,
                                worse_row.hmm_start,
                                worse_row.hmm_end,
                                worse_row.query_start,
                                worse_row.query_end])

        elif len(red_df) == 0:
            print('Empty ORF DF! This should not happen. ORF: ', orf)
            sys.exit()

        else:
            goodrow = red_df.iloc[0]
            if float(goodrow.evalue) < threshold:
                orflist.append([goodrow.family_hmm,
                                genome_id,
                                goodrow.orf_id,
                                goodrow.hmm_length,
                                goodrow.evalue,
                                goodrow.hmm_start,
                                goodrow.hmm_end,
                                goodrow.query_start,
                                goodrow.query_end,
                                'NaN', 'NaN', 'NaN', 'NaN', 'NaN'])


    orf_df = pd.DataFrame(orflist, columns=orflist_header)

    orf_df.to_csv(outdir + '/hmmscan/' + genome_id + '/' + genome_id + '.parse', sep='\t', index=False)

    return outdir + '/hmmscan/' + genome_id + '/' + genome_id + '.parse'

def test():
    args = parser.parse_args()
    args = parser.parse_args()
    nucdir = str(Path(args.nucdir).absolute())
    hmmdir = str(Path(args.hmmdir).absolute())
    outdir = str(args.outdir)
    threshold = float(args.evalue)
    threads = int(args.threads)
    ran_prodigal = args.already_predicted
    already_scanned = args.already_scanned
    no_seqs = args.no_seqs
    best = args.best


    p = Pool(threads)

    # Make output directory
    if not os.path.exists(outdir):
        os.system('mkdir ' + outdir)
        outdir = str(Path(outdir).absolute())
    else:
        outdir = str(Path(outdir).absolute())

    # Get list of paths of all fastas
    fastalist_wpath = list(map(lambda file: os.path.join(nucdir, file), os.listdir(nucdir)))

    # Get list of all fastas
    #fastalist = list(map(lambda file: file.split('.f')[0], os.listdir(fastadir)))

    # Get list of paths of all HMM files
    hmmlist_wpath = list(map(lambda file: os.path.join(hmmdir, file), os.listdir(hmmdir)))

    # Get list of all HMMs
    hmmlist = list(map(lambda file: file.split('.hmm')[0], os.listdir(hmmdir)))

    hmm_outfiles = []




    if not already_scanned:
        if not ran_prodigal:
            #Make folder for proteins
            os.mkdir(outdir + '/proteins')

            #Predict genes for nucleotide fastas
            p.map(lambda x: run_prodigal(x, outdir), fastalist_wpath)

            #Generate binary files for hmmsearch
            hmmpress(hmmlist_wpath, outdir)


        protdir = outdir + '/proteins'

        protlist_wpath = list(map(lambda file: os.path.join(protdir, file), os.listdir(protdir)))

        #Get list of protein files without full path
        protlist = list(map(lambda path: path.split('/')[0], protlist_wpath))

        #Make directory to store hmmsearch outfiles
        if not os.path.exists(outdir + '/hmmscan/'):
            os.system('mkdir ' + outdir + '/hmmscan/')

        def run_hmms(fastafile):
            fastaoutdir = outdir + '/hmmscan/' + fastafile.split('/')[-1].split('.faa')[0].split('.fna')[0].split('.fasta')[0].split('.fa')[0]
            # Make outdir for HMMs
            if not os.path.exists(fastaoutdir):
                os.system('mkdir ' + fastaoutdir)
            #Make symbolic link
            if len(list(filter(lambda x: '.faa' in x, os.listdir(fastaoutdir)))) == 0:
                os.system('ln -s ' + fastafile + ' ' + fastaoutdir + '/')
            hmm_outfiles.append([])

            # Run all HMMs for fastafile
            return run_hmmscan(fastafile, outdir, threshold)

        hmm_outfiles = list(p.map(run_hmms, protlist_wpath))

    #Make sure these variables are loaded in case you activated -already_scanned
    if already_scanned:
        protdir = outdir + '/proteins'

        protlist_wpath = list(map(lambda file: os.path.join(protdir, file), os.listdir(protdir)))

        #Get list of protein files without full path
        protlist = list(map(lambda path: path.split('/')[0], protlist_wpath))

    all_df_list = list(p.map(lambda x: pd.read_csv(x, sep='\t'), hmm_outfiles))
    all_df = pd.concat(all_df_list, sort=False)

    recs_list_by_hmm = extract_hits(all_df, threads, outdir)
    print(all_df.head())






    #Make directory to store fasta hits
    if not os.path.exists(outdir + '/' + 'fastas'):
        os.system('mkdir ' + outdir + '/' + 'fastas')

    make_hitstable_df(recs_list_by_hmm, hmmlist, protlist, outdir)
    print("End make_hitstable_df")
    sys.exit()

    if not no_seqs:
        print("Getting recs and writing to fasta...")

        hmms_written = list(p.map(lambda hits:
                                        write_recs(
                                        #Actual list of recs
                                        hits[0],
                                        #Name of HMM to write (for fasta name)
                                        hits[1],
                                        outdir),
                                        #List of recs to iterate over
                                        recs_list_by_hmm))

    print("Good so far!")

    sys.exit()

    return


def main():
    args = parser.parse_args()

    test()

    if args.nuc:
        nuc_workflow(args)
    elif args.prot:
        prot_workflow(args)
    else:
        print("Please specify a workflow- either nucleotide or protein. Exiting...")
        sys.exit()


if __name__ == "__main__":
    main()
