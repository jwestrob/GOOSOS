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
parser.add_argument('-nucdir', metavar='[NUCLEOTIDE FASTA DIRECTORY]',
                    help='Directory containing nucleotide sequences.')
parser.add_argument('-protdir', metavar='[PRODIGAL PROTEIN FASTA DIRECTORY]',
                    help='Directory of prodigal output for your fastas in amino acid format.', default=None)
parser.add_argument('-hmmdir', metavar='[HMM DIRECTORY]', help="Directory containing HMM files to scan with")
parser.add_argument('-outdir', metavar='[Output Directory]', default='output', help="Directory to store output files")
parser.add_argument('-evalue', metavar='[EVALUE THRESHOLD]', default=0.01,
                    help="Evalue threshold for HMMsearch. Default: 0.01.")
parser.add_argument('-threads', metavar='[NUM THREADS]', default=1, help="Number of threads to use")
parser.add_argument('-already_predicted', default=False, action='store_true', help='For if you already ran Prodigal.')
parser.add_argument('-already_scanned', default=False, action='store_true', help='For if you already ran the HMMs')
parser.add_argument('-no_seqs', default=False, action='store_true', help='Dont pull out sequences to fasta')
parser.add_argument('-best', default=False, action='store_true', help='Only pull out best hit per genome')
parser.add_argument('-align', default=False, action='store_true', help='Align fastas once extracted')
parser.add_argument('-accurate', default=False, action='store_true', help='If aligning, use accurate (SLOW) mafft parameters.')
parser.add_argument('-cut_nc', default=False, action='store_true', help='If using KEGG HMMs, use the --cut_nc option during hmmsearch (built in cutoffs)')
parser.add_argument('-cut_ga', default=False, action='store_true', help='If using PFAM HMMs, use the --cut_ga option during hmmsearch (built in cutoffs)')


def hmmpress(hmmlist_wpath, outdir, cut_nc, cut_ga):
    #Concatenate all hmm files together and press them into a binary
    list_of_hmms = ' '.join(hmmlist_wpath)

    #Make folder to store hmmpress files in
    if not os.path.exists(outdir + '/hmmpress'):
        os.mkdir(outdir + '/hmmpress')

    #Rename the NAME field on each HMM file to be consistent with filename
    #Since hmmscan labels every hit with whatever's in that field
    for hmmfile in hmmlist_wpath:
        rename(hmmfile, outdir + '/hmmpress/')

    if cut_nc or cut_ga:
        hmm_thresh_list = []

        for hmmfile in hmmlist_wpath:
            hmmname = hmmfile.split('.hmm')[0].split('/')[-1]
            with open(hmmfile, 'r') as infile:
                lines = [x.rstrip() for x in infile.readlines()]
            nc = list(filter(lambda x: x.startswith('NC'), lines))

            if len(nc) == 0:
                nc = 0
            else:
                #Grab threshold from middle element
                nc = nc[0].split()[1]

            ga = list(filter(lambda x: x.startswith('GA'), lines))

            if len(ga) == 0:
                ga = 0
            else:
                #Grab threshold from middle element
                ga = ga[0].split()[1]

            thresh  = max(nc, ga)
            hmm_thresh_list.append([hmmname, thresh])

        hmm_thresh_dict = dict(hmm_thresh_list)
    else:
        hmm_thresh_list = None
        hmm_thresh_dict = None




    hmmpressdir = outdir + '/hmmpress/'

    print('cat ' + hmmpressdir + '*.hmm > ' + hmmpressdir + 'concatenated_hmms.hmm')

    os.system('cat ' + hmmpressdir + '*.hmm > ' + hmmpressdir + 'concatenated_hmms.hmm')

    os.system('hmmpress ' + hmmpressdir + 'concatenated_hmms.hmm')

    return hmm_thresh_dict

def rename(hmmfile, hmmdir):
    with open(hmmfile, 'r') as f:
        lines = f.readlines()
    name_line = lines[1].split()
    name_line[1] = hmmfile.split('.hmm')[0].split('/')[-1]
    lines[1] = ('\t').join(name_line)
    with open(hmmdir + hmmfile.split('/')[-1], 'w') as f:
        for item in lines:
            f.write("%s\n" % item)
    return

def run_hmmscan(protfile, outdir, threshold, best, cut_nc, cut_ga):
    genome_id = protfile.split('/')[-1].split('.faa')[0].split('.fna')[0].split('.fa')[0].split('.fasta')[0]

    if len(list(filter(lambda x: '_hmmsearch.out' in x, os.listdir(outdir + '/hmmscan/' + '/')))) > 0:
        return parse_hmmdomtbl(outdir, genome_id + '_hmmsearch.out', threshold, best)

    #print(protein_id, hmmfile)
    if not cut_nc and not cut_ga:
        cmd = 'hmmscan --domtblout ' + outdir + '/hmmscan/' + genome_id + '/' + genome_id + '_hmmsearch.out  --notextw --cpu ' \
                + str(1) + ' ' + outdir + '/hmmpress/concatenated_hmms.hmm ' + protfile + ' > /dev/null 2>&1'
    elif cut_nc:
        cmd = 'hmmscan --domtblout ' + outdir + '/hmmscan/' + genome_id + '/' + genome_id + '_hmmsearch.out  --notextw --cut_nc --cpu ' \
                + str(1) + ' ' + outdir + '/hmmpress/concatenated_hmms.hmm ' + protfile + ' > /dev/null 2>&1'
    elif cut_ga:
        cmd = 'hmmscan --domtblout ' + outdir + '/hmmscan/' + genome_id + '/' + genome_id + '_hmmsearch.out  --notextw --cut_ga --cpu ' \
                + str(1) + ' ' + outdir + '/hmmpress/concatenated_hmms.hmm ' + protfile + ' > /dev/null 2>&1'
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
    return parse_hmmdomtbl(outdir, genome_id + '_hmmsearch.out', threshold, best)


def get_rec_for_hit(genome_id, orf, outdir):
    genome_id = str(genome_id)
    genome_dir = outdir + '/hmmscan/' + genome_id + '/'
    protfile = list(filter(lambda x: '.faa' in x, os.listdir(genome_dir)))[0]
    genome_recs = list(SeqIO.parse(genome_dir + protfile, 'fasta'))

    desired_hit = list(filter(lambda x: orf == x.id, genome_recs))[0]
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

    #Make sure that you extract the best hit
    #dedupe_df = dedupe(all_df)

    #I could use a map, but like... why
    for hmm in all_df['family_hmm'].unique().tolist():
        red_df_nothresh = all_df[all_df['family_hmm'] == hmm]
        red_df = red_df_nothresh[red_df_nothresh['above_threshold']]
        recs_by_hmm.append([extract_hits_by_hmm(red_df, threads, outdir), hmm])

    return recs_by_hmm

def make_hitstable_df(recs_by_hmm, hmmlist, fastalist, outdir):

    # Make matrix of zeros to store hits
    print("Making hits matrix...")
    hitstable = np.zeros((len(hmmlist), len(fastalist)))

    hits = pd.DataFrame(hitstable).T
    hits.columns = hmmlist

    hits.index = fastalist

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

            hits[hmm][genome_id] += 1


    #Write it to tsv in outdir without index (annoying)
    colnames = hits.columns.values.tolist()

    hits.to_csv(outdir + '/HITSTABLE.tsv', sep='\t')

def write_recs(recs_for_hmm, hmm_name, outdir):
    fasta_outdir = outdir + '/fastas'

    recs_for_hmm = list(filter(lambda x: type(x) is not None, recs_for_hmm))

    print("Writing recs for " + hmm_name)


    SeqIO.write(recs_for_hmm, fasta_outdir + '/' + hmm_name + '_hits.faa', 'fasta')
    return

def run_prodigal(fastafile_wpath, outdir):
    #print('prodigal -i '+ fastafile_wpath + ' -a ' + outdir + '/proteins/' +
    #                        fastafile_wpath.split('/')[-1] + '.faa -m -p single > /dev/null 2>&1')
    os.system('prodigal -i '+ fastafile_wpath + ' -a ' + outdir + '/proteins/' +
                            fastafile_wpath.split('/')[-1] + '.faa -m -p single > /dev/null 2>&1')

    protdir = outdir + '/proteins/'
    if os.path.getsize(protdir + fastafile_wpath.split('/')[-1] + '.faa') == 0:
        print('No genes predicted for ' + fastafile_wpath.split('/')[-1])
        os.system('rm ' + protdir + fastafile_wpath.split('/')[-1] + '.faa')
    else:
        print('Genes predicted for ' + fastafile_wpath.split('/')[-1])
    return

def parse_hmmdomtbl(outdir, hmmoutfile, threshold, best):

    """
    Takes output file from HMMscan (--domtblout), parses it, and yields a
    genome-specific dataframe of each hit (above the evalue threshold) for each
    domain.
    Now contains bitscore information as well.
    """

    genome_id = hmmoutfile.split('_hmmsearch.out')[0].split('.fasta')[0].split('.fna')[0].split('.fa')[0]
    hmmoutfile_wpath = outdir + '/hmmscan/' + genome_id + '/' + hmmoutfile

    with open(hmmoutfile_wpath, 'r') as infile:
        lines = infile.readlines()

    domtbl_header = ['target name', 'accession', 'tlen', 'query name', 'accession', 'qlen', 'E-value',
                    'full_score',  'full_bias',   'dom_#',  'dom_of',  'c-Evalue',  'i-Evalue',  'dom_score',  'dom_bias',  'hmm_from',
                    'hmm_to',  'ali_from', 'ali_to',  'env_from', 'env_to',  'mean_posterior', 'description of target']

    desired_header = ['family_hmm', 'hmm_length', 'orf_id',
                      'query_length', 'bitscore', 'evalue', 'hmm_start',
                      'hmm_end', 'query_start', 'query_end']

    #Remove all lines with '#' beginning character
    lines_filtered = list(filter(lambda x: x[0] != '#', lines))
    if len(lines_filtered) == 0:
        print("No hits for " + genome_id)
        orflist_header = ['family_hmm', 'hmm_length', 'orf_id', 'overall_bitscore', 'overall_evalue', 'dom1_evalue', 'dom1_hmmstart',
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

    #Make pandas DF to store lines, then add column names
    try:
	    lines_df = pd.DataFrame(lines_filtered, columns=domtbl_header)
    except:
           print("Error parsing hmmdomtbl for: ", genome_id)
           sys.exit()
    #Make DF to store properly arranged data
    goodheader_df = pd.DataFrame(columns=desired_header)

    goodheader_df['family_hmm'] = lines_df['target name']
    #Insert genome ID to avoid confusion
    goodheader_df['genome_id'] = pd.Series([genome_id]*len(lines_df))
    goodheader_df['hmm_length'] = lines_df['tlen']
    goodheader_df['orf_id'] = lines_df['query name']
    goodheader_df['query_length'] = lines_df['qlen']
    goodheader_df['bitscore'] = lines_df['full_score']
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
    orflist_header = ['family_hmm', 'genome_id', 'orf_id', 'hmm_length', 'query_length', 'overall_bitscore', 'overall_evalue', 'dom1_cevalue', 'dom1_hmmstart',
                      'dom1_hmmend', 'dom1_querystart', 'dom1_queryend', 'dom2_cevalue', 'dom2_hmmstart',
                      'dom2_hmmend', 'dom2_querystart', 'dom2_queryend']
    for orf in unique_orfs:
        red_df = goodheader_df[goodheader_df['orf_id'] == orf]

        #For the unlikely scenario where you have more than one HMM hit on a single ORF xfor ONE HMM
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
                                goodrow.query_length,
                                goodrow.bitscore,
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
                                goodrow.query_length,
                                goodrow.bitscore,
                                goodrow.evalue,
                                goodrow.c_evalue,
                                goodrow.hmm_start,
                                goodrow.hmm_end,
                                goodrow.query_start,
                                goodrow.query_end,
                                'NaN', 'NaN', 'NaN', 'NaN', 'NaN'])


    orf_df = pd.DataFrame(orflist, columns=orflist_header)
    orf_df = orf_df[orf_df['overall_evalue'].astype(float) <= threshold]
    orf_df = orf_df[orf_df['query_length'].astype(float) >= 0.75*orf_df['hmm_length'].astype(float)]
    #and float(goodheader_df.query_length) > 0.75*float(goodheader_df.hmm_length)

    #Great. Now let's deduplicate.
    if best:
        hmms = orf_df.family_hmm.unique()
        for hmm in hmms:
            red_df = orf_df[orf_df['family_hmm'] == hmm]
            if len(red_df) == 1:
                continue
            else:
                red_df.sort_values(by='overall_bitscore')
                to_drop = []
                for i in range(1, len(red_df)):
                    to_drop.append(red_df.iloc[i].name)
                orf_df.drop(to_drop, axis=0)



    orf_df.to_csv(outdir + '/hmmscan/' + genome_id + '/' + genome_id + '.parse', sep='\t', index=False)

    return outdir + '/hmmscan/' + genome_id + '/' + genome_id + '.parse'

def align_fn(fastafile_wpath, outdir, threads, accurate):
    fastafile_id = fastafile_wpath.split('/')[-1].split('.faa')[0]
    if accurate:
        #print('mafft --localpair --thread ' + str(threads) + ' --maxiterate 1000 ' + fastafile_wpath + ' > ' + outdir + '/alignments/' + fastafile_id + '_ALN.mfaa')
        os.system('mafft --localpair --thread ' + str(threads) + ' --maxiterate 1000 ' + fastafile_wpath + ' > '
                + outdir + '/alignments/' + fastafile_id + '_ALN.mfaa')
    else:
        os.system('mafft --thread ' + str(threads) + ' ' + fastafile_wpath + ' > '
                + outdir + '/' + alignments + '/' + fastafile_id + '_ALN.mfaa')
    return


def fetch_outfiles(outdir, threshold, threads, best):
    p = Pool(threads)
    hmmscandir = os.path.join(outdir, 'hmmscan')
    subdirectories = os.listdir(hmmscandir)
    subdirectories_wpath = list(map(lambda dir: os.path.join(hmmscandir, dir), subdirectories))

    def get_outfile(dir):
        all_files = os.listdir(dir)
        hmmoutfile = list(filter(lambda x: '_hmmsearch.out' in x, all_files))[0]
        return hmmoutfile

    hmmoutfiles = list(map(get_outfile, subdirectories_wpath))

    parsed_hmm_outfiles = list(p.map(lambda outfile: parse_hmmdomtbl(outdir, outfile, threshold, best),
                                     hmmoutfiles))
    return parsed_hmm_outfiles


def run_hmms(fastafile, outdir, threshold, best, cut_nc, cut_ga):

    fastaoutdir = outdir + '/hmmscan/' + fastafile.split('/')[-1].split('.faa')[0].split('.fna')[0].split('.fasta')[0].split('.fa')[0]

    # Make outdir for HMMs
    if not os.path.exists(fastaoutdir):
        os.system('mkdir ' + fastaoutdir)
    #Make symbolic link
    if len(list(filter(lambda x: '.faa' in x, os.listdir(fastaoutdir)))) == 0:
        os.system('ln -s ' + fastafile + ' ' + fastaoutdir + '/')

    # Run all HMMs for fastafile
    return run_hmmscan(fastafile, outdir, threshold, best, cut_nc, cut_ga)

def mark_with_threshold(all_df, hmm_thresh_dict):

    if hmm_thresh_dict is not  None:
        all_df['above_threshold'] = all_df.apply(lambda x: x.overall_bitscore >= float(hmm_thresh_dict[x.family_hmm]), axis=1)

    return all_df

def main():
    args = parser.parse_args()
    if args.nucdir is not None:
        nucdir = str(Path(args.nucdir).absolute())
    else:
        nucdir = None
    if args.protdir is not None:
        prodigaldir = str(Path(args.protdir).absolute())
    else:
        prodigaldir = None

    if nucdir is None and prodigaldir is None:
        print("You need to provide a protein or nucleotide directory.")
        sys.exit(420)
    hmmdir = str(Path(args.hmmdir).absolute())
    outdir = str(args.outdir)
    threshold = float(args.evalue)
    threads = int(args.threads)
    ran_prodigal = args.already_predicted
    already_scanned = args.already_scanned
    no_seqs = args.no_seqs
    best = args.best
    align = args.align
    accurate = args.accurate
    cut_nc = args.cut_nc
    cut_ga = args.cut_ga


    have_proteins = True if prodigaldir is not None else False

    p = Pool(threads)

    # Make output directory
    if not os.path.exists(outdir):
        os.system('mkdir ' + outdir)
        outdir = str(Path(outdir).absolute())
    else:
        outdir = str(Path(outdir).absolute())

    if not have_proteins:
        # Get list of paths of all fastas
        fastalist_wpath = list(map(lambda file: os.path.join(nucdir, file), os.listdir(nucdir)))

    # Get list of all fastas
    #fastalist = list(map(lambda file: file.split('.f')[0], os.listdir(fastadir)))

    # Get list of paths of all HMM files
    hmmlist_wpath = list(map(lambda file: os.path.join(hmmdir, file), os.listdir(hmmdir)))

    # Get list of all HMMs
    hmmlist = list(map(lambda file: file.split('.hmm')[0], os.listdir(hmmdir)))

    if not already_scanned:
        #Generate binary files for hmmsearch
        hmm_thresh_dict = hmmpress(hmmlist_wpath, outdir, cut_nc, cut_ga)
        np.save(os.path.join(outdir, 'hmm_thresh_dict.npy'), hmm_thresh_dict)
    else:
        np.load(os.path.join(outdir, 'hmm_thresh_dict.npy'), allow_pickle=True)

    parsed_hmm_outfiles = []


    if not already_scanned:
        if not ran_prodigal and not have_proteins:
            #Make folder for proteins
            if not os.path.exists(outdir + '/proteins'):
                os.mkdir(outdir + '/proteins')

            #Predict genes for nucleotide fastas
            p.map(lambda x: run_prodigal(x, outdir), fastalist_wpath)


        if not have_proteins:
            protdir = outdir + '/proteins'
        else:
            protdir = prodigaldir

        protlist_wpath = list(map(lambda file: os.path.join(protdir, file), os.listdir(protdir)))

        #Get list of protein files without full path
        protlist = list(map(lambda path: path.split('/')[-1].split('.fna')[0].split('.fa')[0].split('.fasta')[0],
                        protlist_wpath))

        #Make directory to store hmmsearch outfiles
        if not os.path.exists(outdir + '/hmmscan/'):
            os.system('mkdir ' + outdir + '/hmmscan/')

        if cut_nc is False and cut_ga is False:
            print('hmmscan will be run with no built-in thresholds. Please make sure you filter afterwards if you want quality filtering.')
            print('All you need to do is delete unsuitable rows from all_hits_evalues_df and rerun with the -already_scanned flag.')

        #Make sure you get rid of any Nones
        parsed_hmm_outfiles = list(filter(lambda x: x is not None, list(p.map(lambda x: run_hmms(x, outdir, threshold, best, cut_nc, cut_ga), protlist_wpath))))

        all_df_list = list(p.map(lambda x: pd.read_csv(x, sep='\t'), parsed_hmm_outfiles))
        all_df_init = pd.concat(all_df_list, sort=False)

        all_df_thresh = mark_with_threshold(all_df_init, hmm_thresh_dict)

        all_df = all_df_thresh[all_df_thresh['above_threshold']]

        all_df.to_csv(outdir + '/all_hits_evalues_df.tsv', sep='\t', index=False)
        print("Hits information written to all_hits_evalues_df.tsv.")
    #Make sure these variables are loaded in case you activated -already_scanned
    if already_scanned:
        if not have_proteins:
            protdir = outdir + '/proteins'
        else:
            protdir = prodigaldir

        if not os.path.exists(outdir + '/all_hits_evalues_df.tsv'):
            protlist_wpath = list(map(lambda file: os.path.join(protdir, file), os.listdir(protdir)))

            #Make sure you get rid of any Nones
            parsed_hmm_outfiles = list(filter(lambda x: x is not None, list(p.map(lambda x: run_hmms(x, outdir, threshold, best, cut_nc, cut_ga), protlist_wpath))))

            all_df_list = list(p.map(lambda x: pd.read_csv(x, sep='\t'), parsed_hmm_outfiles))
            all_df = pd.concat(all_df_list, sort=False)

            all_df.to_csv(outdir + '/all_hits_evalues_df.tsv', sep='\t', index=False)

        protlist_wpath = list(map(lambda file: os.path.join(protdir, file), os.listdir(protdir)))

        #Get list of protein files without full path
        protlist = list(map(lambda path: path.split('/')[-1].split('.fna')[0].split('.fa')[0].split('.fasta')[0],
                        protlist_wpath))

        all_df = pd.read_csv(outdir + '/all_hits_evalues_df.tsv', sep='\t')


    recs_list_by_hmm = extract_hits(all_df, threads, outdir)

    #Make directory to store fasta hits
    if not os.path.exists(outdir + '/' + 'fastas'):
        os.system('mkdir ' + outdir + '/' + 'fastas')

    make_hitstable_df(recs_list_by_hmm, hmmlist, protlist, outdir)

    if not no_seqs:
        print("Getting recs and writing to fasta...")

        hmms_written = list(map(lambda hits:
                                        write_recs(
                                        #Actual list of recs
                                        hits[0],
                                        #Name of HMM to write (for fasta name)
                                        hits[1],
                                        outdir),
                                        #List of recs to iterate over
                                        recs_list_by_hmm))

        def sort_fasta(fastafile_wpath, fastalist):
            in_recs = list(SeqIO.parse(fastafile_wpath, 'fasta'))
            out_recs = []
            for fasta in fastalist:
                rec_to_append = list(filter(lambda x: fasta in x.id, in_recs))
                out_recs.append(rec_to_append)
            SeqIO.write(out_recs, fastafile_wpath, 'fasta')
            return

        if align:
            out_fastas = os.listdir(outdir + '/fastas')
            out_fastas = list(map(lambda x: os.path.join(outdir + '/fastas', x), out_fastas))
            #Sort those fastas before aligning!
            p.map(lambda x: sort_fasta(x, protlist), out_fastas)
            os.system('mkdir ' + outdir + '/alignments')
            list(map(lambda x: align_fn(x, outdir, threads, accurate), out_fastas))


    prot_series = pd.Series(protlist)
    prot_series.to_csv(outdir + '/genome_order.txt', sep=',', index=False, header=None)

    print("Order of sorted fastas written to genome_order.txt.")

    print("You did it!")

if __name__ == "__main__":
    main()
