from pathos.multiprocessing import ProcessingPool as Pool
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
#from Bio.Alphabet import IUPAC
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
parser.add_argument('-synteny', default=False, action='store_true', help='If doing RP16, ensure syntenic blocks (>=8 hits on single contigs per genome)')

#------------------------
#      FUNCTION ZOO
#------------------------

def gather_hmms(hmmlist_wpath, outdir, cut_nc, cut_ga):
    #Concatenate all hmm files together and press them into a binary
    list_of_hmms = ' '.join(hmmlist_wpath)

    #Make folder to store hmm files in
    if not os.path.exists(outdir + '/renamed_hmms'):
        os.mkdir(outdir + '/renamed_hmms')

    #Rename the NAME field on each HMM file to be consistent with filename
    #Since hmmsearch labels every hit with whatever's in that field
    for hmmfile in hmmlist_wpath:
        rename(hmmfile, outdir + '/renamed_hmms/')

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
                nc = float(nc[0].split()[1])

            ga = list(filter(lambda x: x.startswith('GA'), lines))

            if len(ga) == 0:
                ga = 0
            else:
                #Grab threshold from middle element
                ga = float(ga[0].split()[1])

            thresh  = max(float(nc), float(ga))
            hmm_thresh_list.append([hmmname, thresh])

        hmm_thresh_dict = dict(hmm_thresh_list)
    else:
        hmm_thresh_dict = None

    renamed_dir = outdir + '/renamed_hmms/'

    print('cat ' + renamed_dir + '*.hmm > ' + renamed_dir + 'concatenated_hmms.hmm')

    os.system('cat ' + renamed_dir + '*.hmm > ' + renamed_dir + 'concatenated_hmms.hmm')

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

def run_hmmsearch(protfile, outdir, threshold, best, cut_nc, cut_ga):
    genome_id = protfile.split('/')[-1].split('.faa')[0].split('.fna')[0].split('.fa')[0].split('.fasta')[0]
    #if len(list(filter(lambda x: '_hmmsearch.out' in x, os.listdir(outdir + '/hmmsearch/' + '/')))) > 0:
    #    return parse_hmmdomtbl(outdir, genome_id + '_hmmsearch.out', threshold, best)

    #print(protein_id, hmmfile)
    if not cut_nc and not cut_ga:
        cmd = 'hmmsearch --domtblout ' + outdir + '/hmmsearch/' + genome_id + '/' + genome_id + '_hmmsearch.out  --notextw --cpu ' \
                + str(1) + ' -E ' + str(threshold) + ' ' + outdir + '/renamed_hmms/concatenated_hmms.hmm ' + protfile + ' > /dev/null 2>&1'
    elif cut_nc:
        cmd = 'hmmsearch --domtblout ' + outdir + '/hmmsearch/' + genome_id + '/' + genome_id + '_hmmsearch.out  --notextw --cut_nc --cpu ' \
                + str(1) + ' ' + outdir + '/renamed_hmms/concatenated_hmms.hmm ' + protfile + ' > /dev/null 2>&1'
    elif cut_ga:
        cmd = 'hmmsearch --domtblout ' + outdir + '/hmmsearch/' + genome_id + '/' + genome_id + '_hmmsearch.out  --notextw --cut_ga --cpu ' \
                + str(1) + ' ' + outdir + '/renamed_hmms/concatenated_hmms.hmm ' + protfile + ' > /dev/null 2>&1'
    #print(cmd)
    result = subprocess.run(cmd, shell=True, check=True)
    if result.returncode != 0:
        print(result)
        print('HMMsearch error (check for empty sequences in your protein FASTAs)')
        print('genome_id: ', genome_id)
        #print('hmmfile: ', hmmfile)
        sys.exit()
    #print(result)
    #Parse file with awk/perl nonsense; generate .parse file
    parse_hmmdomtbl_multidomain(outdir, genome_id + '_hmmsearch.out', threshold)
    return parse_hmmdomtbl(outdir, genome_id + '_hmmsearch.out', threshold, best)

def get_rec_for_hit(genome_id, orf, outdir):
    genome_id = str(genome_id)
    genome_dir = outdir + '/hmmsearch/' + genome_id + '/'
    protfile = list(filter(lambda x: '.faa' in x or '.fa' in x, os.listdir(genome_dir)))[0]
    genome_generator = SeqIO.parse(genome_dir + protfile, 'fasta')
    try:
        desired_hit = list(filter(lambda x: orf == x.id, genome_generator))[0]
    except IndexError:
        return None
    if not desired_hit.id.startswith(genome_id + '|'):
        desired_hit.id = genome_id + '|' + desired_hit.id

    return desired_hit

def extract_hits_by_hmm(red_df, threads, outdir):
    print("Extracting hits for " + red_df.iloc[0].family_hmm)
    p2 = Pool(threads)

    #Make list of genome_id / orf pairs
    id_orf_list = red_df[['genome_id', 'orf_id']].values.tolist()

    recs = list(p2.map(lambda id_orf: get_rec_for_hit(id_orf[0], id_orf[1], outdir), id_orf_list))
    recs = list(filter(lambda x: x is not None, recs))
    return recs

def extract_hits(all_df, threads, outdir):
    #List of recs (value to return)
    recs_by_hmm = []

    #Make sure that you extract the best hit
    #dedupe_df = dedupe(all_df)

    #I could use a map, but like... why
    for hmm in all_df['family_hmm'].unique().tolist():
        if 'above_threshold' in all_df.columns.tolist():
            red_df_nothresh = all_df[all_df['family_hmm'] == hmm]
            red_df = red_df_nothresh[red_df_nothresh['above_threshold']]
        else:
            red_df = all_df[all_df['family_hmm'] == hmm]
        recs_by_hmm.append([extract_hits_by_hmm(red_df, threads, outdir), hmm])

    return recs_by_hmm

def extract_hits_2(all_df, threads, outdir):
    #Takes all_df and generates recs for each HMM by extracting on a
    #file-by-file basis (bc of disk I/O limitations with large files)

    recs_by_file = []
    recs_by_hmm = []


    def grab_recs_by_genome(genome_id, outdir=outdir):
        #Get all ORFs corresponding to that genome/input file
        desired_orfs = all_df[all_df.genome_id == genome_id].orf_id.tolist()
        genome_dir = outdir + '/hmmsearch/' + genome_id + '/'
        #Grab protein file from directory within hmmsearch folder
        protfile = list(filter(lambda x: x.endswith('.faa') or x.endswith('.fa'), os.listdir(genome_dir)))[0]

        genome_recs = list(filter(lambda x: x.id in desired_orfs, SeqIO.parse(os.path.join(genome_dir, protfile), 'fasta')))
        return genome_recs

    p2 = Pool(threads)
    recs_by_file = list(p2.map(grab_recs_by_genome,  all_df.genome_id.unique().tolist()))


    flatten = lambda l: [item for sublist in l for item in sublist]

    all_recs = flatten(recs_by_file)

    for hmm in all_df.family_hmm.unique().tolist():
        desired_orfs_2 = all_df[all_df.family_hmm == hmm].orf_id.tolist()
        hmm_recs = list(filter(lambda x: x.id in desired_orfs_2, all_recs))
        recs_by_hmm.append([hmm_recs, hmm])

    return recs_by_hmm

def extract_hits_3(all_df, threads, outdir):
    os.system('mkdir ' + os.path.join(outdir, 'pullseq_tmp'))

    #Takes all_df and generates recs for each HMM by extracting with pullseq

    genome_tmp = os.path.join(outdir, 'pullseq_tmp/genome_fastas')
    os.system('mkdir ' + genome_tmp)
    orfids_tmp = os.path.join(outdir, 'pullseq_tmp/orfids_by_genome')
    os.system('mkdir ' + orfids_tmp)


    recs_by_file = []
    recs_by_hmm = []

    def grab_recs_by_genome(genome_id, outdir=outdir, genome_tmp=genome_tmp, orfids_tmp=orfids_tmp):
        #Get all ORFs corresponding to that genome/input file
        desired_orfs = all_df[all_df.genome_id == genome_id].orf_id.tolist()
        genome_dir = outdir + '/hmmsearch/' + genome_id + '/'
        #Grab protein file from directory within hmmsearch folder
        protfile = list(filter(lambda x: x.endswith('.faa') or x.endswith('.fa'), os.listdir(genome_dir)))[0]

        with open(os.path.join(orfids_tmp, genome_id + '.txt'), 'w') as outfile:
            for element in desired_orfs:
                outfile.writelines(str(element) + '\n')

        os.system('pullseq -i ' + os.path.join(genome_dir, protfile) + ' -n ' + os.path.join(orfids_tmp, genome_id + '.txt > ' + genome_tmp + '/' + genome_id + 'allhits.faa'))
        genome_recs = list(SeqIO.parse(os.path.join(genome_tmp, genome_id + 'allhits.faa'), 'fasta'))
        return genome_recs

    p2 = Pool(threads)
    print("Fetching orfs...")
    recs_by_file = list(p2.map(grab_recs_by_genome,  all_df.genome_id.unique().tolist()))


    flatten = lambda l: [item for sublist in l for item in sublist]

    all_recs = pd.Series(flatten(recs_by_file))
    print("Separating orfs...")
    def separate_orfs(orflist, hmmname, all_df=all_df, all_recs=all_recs):
        #all_recs is pd.Series now
        hmm_recs = all_recs[all_recs.apply(lambda x: x.id in orflist)]
        return([hmm_recs, hmmname])

    #This is such a dope line of code. Swag
    recs_by_hmm = list(p2.map(
                    lambda hmm: separate_orfs(all_df[all_df.family_hmm == hmm].orf_id.tolist(), hmm),
                    all_df.family_hmm.tolist()))

    return recs_by_hmm

def extract_hits_4(all_df, threads, protdir, outdir):
    orfids_tmp = os.path.join(outdir, 'pullseq_tmp')
    os.system('mkdir ' + orfids_tmp)
    os.system('mkdir ' + os.path.join(outdir, 'fastas'))
    #Takes all_df and generates recs for each HMM by extracting with pullseq


    def grab_recs_by_hmm(hmm, protdir=protdir, orfids_tmp=orfids_tmp, all_df=all_df):
        hits_fasta = os.path.join(outdir + '/fastas',  hmm + '_hits.faa')
        #Create/overwrite existing fasta file for hits
        os.system('>' + hits_fasta)
        idfile = os.path.join(orfids_tmp, hmm + '.txt')
        with open(idfile, 'w') as outfile:
            for orfid in all_df[all_df.family_hmm == hmm].orf_id.tolist():
                outfile.write(str(orfid) + '\n')

        for genome in all_df.genome_id.unique().tolist():

            genome_file = list(filter(lambda x: x.split('.fa')[0].split('.fna')[0] == genome, os.listdir(protdir)))[0]
            genome_file = os.path.join(protdir, genome_file)
            pullseq_cmd = 'cat ' + idfile + ' | pullseq -i ' + genome_file + ' -N >> ' + hits_fasta
            os.system(pullseq_cmd)

        return [all_df[all_df.family_hmm == hmm].orf_id.tolist(), hmm]

    p2 = Pool(threads)
    print("Fetching orfs...")
    rec_ids_by_hmm = list(p2.map(grab_recs_by_hmm,  all_df.family_hmm.unique().tolist()))
    os.system('rm -rf ' + orfids_tmp)


    return rec_ids_by_hmm


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
    Takes output file from HMMsearch (--domtblout), parses it, and yields a
    genome-specific dataframe of each hit (above the evalue threshold) for each
    domain.
    Now contains bitscore information as well.
    """

    genome_id = hmmoutfile.split('_hmmsearch.out')[0].split('.fasta')[0].split('.fna')[0].split('.fa')[0]
    hmmoutfile_wpath = outdir + '/hmmsearch/' + genome_id + '/' + hmmoutfile

    with open(hmmoutfile_wpath, 'r') as infile:
        lines = infile.readlines()

    domtbl_header = ['orf_id', 'accession', 'seqlen', 'family_hmm',
                'accession', 'hmm_len', 'overall_evalue', 'overall_bitscore', 'bias',
                'num_domains', 'num_domains_2', 'dom_cond_evalue', 'dom_ind_evalue',
                'dom_bitscore', 'dom_bias', 'hmm_start', 'hmm_end', 'seq_start', 'seq_end',
                'env_start', 'env_end', 'acc']

    desired_header = ['family_hmm', 'hmm_length', 'orf_id',
                      'query_length', 'bitscore', 'evalue', 'num_domains', 'hmm_start',
                      'hmm_end', 'query_start', 'query_end']

    #Remove all lines with '#' beginning character
    lines_filtered = list(filter(lambda x: x[0] != '#', lines))
    if len(lines_filtered) == 0:
        print("No hits for " + genome_id)
        empty_df = pd.DataFrame(columns = desired_header)
        empty_df.to_csv(outdir + '/hmmsearch/' + genome_id + '/' + genome_id + '.parse', sep='\t', index=False)
        return outdir + '/hmmsearch/' + genome_id + '/' + genome_id + '.parse'

    #Remove newline characters and split by whitespace
    lines_filtered = list(map(lambda x: x.strip('\n').split(), lines_filtered))

    #Python, by default, splits the description at the end, causing dimension mismatch.
    #Let's get rid of the extra elements.
    lines_filtered = list(map(lambda x: x[0:22], lines_filtered))

    #Make pandas DF to store lines, then add column names
    try:
	    lines_df = pd.DataFrame(lines_filtered, columns=domtbl_header)
    except:
           print("Error parsing hmmdomtbl for: ", genome_id)
           sys.exit()
    #Make DF to store properly arranged data
    goodheader_df = pd.DataFrame(columns=desired_header)

    goodheader_df['orf_id'] = lines_df['orf_id']
    #Insert genome ID to avoid confusion
    goodheader_df['genome_id'] = pd.Series([genome_id]*len(lines_df))
    goodheader_df['hmm_length'] = lines_df['hmm_len']
    goodheader_df['num_domains'] = lines_df['num_domains_2']
    goodheader_df['family_hmm'] = lines_df['family_hmm']
    goodheader_df['query_length'] = lines_df['seqlen']
    goodheader_df['bitscore'] = lines_df['overall_bitscore']
    goodheader_df['evalue'] = lines_df['overall_evalue']
    goodheader_df['c_evalue'] = lines_df['dom_cond_evalue']
    goodheader_df['hmm_start'] = lines_df['hmm_start']
    goodheader_df['hmm_end'] = lines_df['hmm_end']
    goodheader_df['query_start'] = lines_df['seq_start']
    goodheader_df['query_end'] = lines_df['seq_end']

    unique_orfs = goodheader_df['orf_id'].unique()
    #print("Unique orfs: ")
    #print(unique_orfs)

    orflist = []
    orflist_header = ['family_hmm', 'genome_id', 'orf_id', 'hmm_length', 'query_length', 'num_domains', 'overall_bitscore', 'overall_evalue', 'dom1_cevalue', 'dom1_hmmstart',
                      'dom1_hmmend', 'dom1_querystart', 'dom1_queryend', 'dom2_cevalue', 'dom2_hmmstart',
                      'dom2_hmmend', 'dom2_querystart', 'dom2_queryend']
    for orf in unique_orfs:
        red_df = goodheader_df[goodheader_df['orf_id'] == orf]

        #For the unlikely scenario where you have more than one HMM hit on a single ORF for ONE HMM
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
                                goodrow.num_domains,
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
                                goodrow.num_domains,
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



    orf_df.to_csv(outdir + '/hmmsearch/' + genome_id + '/' + genome_id + '.parse', sep='\t', index=False)

    return outdir + '/hmmsearch/' + genome_id + '/' + genome_id + '.parse'

def parse_hmmdomtbl_multidomain(outdir, hmmoutfile, threshold):

    """
    Takes output file from HMMsearch (--domtblout), parses it, and yields a
    genome-specific dataframe of each hit (above the evalue threshold) for each
    domain.
    """

    genome_id = hmmoutfile.split('_hmmsearch.out')[0].split('.fasta')[0].split('.fna')[0].split('.fa')[0]
    hmmoutfile_wpath = outdir + '/hmmsearch/' + genome_id + '/' + hmmoutfile

    with open(hmmoutfile_wpath, 'r') as infile:
        lines = infile.readlines()

    domtbl_header = ['orf_id', 'accession', 'seqlen', 'family_hmm',
                'accession', 'hmm_len', 'overall_evalue', 'overall_bitscore', 'bias',
                'num_domains', 'num_domains_2', 'dom_cond_evalue', 'dom_ind_evalue',
                'dom_bitscore', 'dom_bias', 'hmm_start', 'hmm_end', 'seq_start', 'seq_end',
                'env_start', 'env_end', 'acc']

    desired_header = ['family_hmm', 'hmm_length', 'orf_id',
                      'query_length', 'bitscore', 'evalue', 'num_domains', 'hmm_start',
                      'hmm_end', 'query_start', 'query_end']

    #Remove all lines with '#' beginning character
    lines_filtered = list(filter(lambda x: x[0] != '#', lines))
    if len(lines_filtered) == 0:
        print("No hits for " + genome_id)
        empty_df = pd.DataFrame(columns = desired_header)
        empty_df.to_csv(outdir + '/hmmsearch/' + genome_id + '/' + genome_id + '.parse', sep='\t', index=False)
        return outdir + '/hmmsearch/' + genome_id + '/' + genome_id + '.parse'

    #Remove newline characters and split by whitespace
    lines_filtered = list(map(lambda x: x.strip('\n').split(), lines_filtered))

    #Python, by default, splits the description at the end, causing dimension mismatch.
    #Let's get rid of the extra elements.
    lines_filtered = list(map(lambda x: x[0:22], lines_filtered))

    #Make pandas DF to store lines, then add column names
    try:
	    lines_df = pd.DataFrame(lines_filtered, columns=domtbl_header)
    except:
           print("Error parsing hmmdomtbl for: ", genome_id)
           sys.exit()
    #Make DF to store properly arranged data
    goodheader_df = pd.DataFrame(columns=desired_header)

    goodheader_df['orf_id'] = lines_df['orf_id']
    #Insert genome ID to avoid confusion
    goodheader_df['genome_id'] = pd.Series([genome_id]*len(lines_df))

    goodheader_df['hmm_length'] = lines_df['hmm_len']
    goodheader_df['num_domains'] = lines_df['num_domains_2']
    goodheader_df['family_hmm'] = lines_df['family_hmm']
    goodheader_df['query_length'] = lines_df['seqlen']
    goodheader_df['bitscore'] = lines_df['overall_bitscore']
    goodheader_df['overall_evalue'] = lines_df['overall_evalue']
    goodheader_df['dom_c_evalue'] = lines_df['dom_cond_evalue']
    goodheader_df['hmm_start'] = lines_df['hmm_start']
    goodheader_df['hmm_end'] = lines_df['hmm_end']
    goodheader_df['query_start'] = lines_df['seq_start']
    goodheader_df['query_end'] = lines_df['seq_end']

    goodheader_df = goodheader_df[goodheader_df['overall_evalue'].astype(float) <= threshold]
    goodheader_df = goodheader_df[goodheader_df['query_length'].astype(float) >= 0.75*goodheader_df['hmm_length'].astype(float)]

    goodheader_df.to_csv(outdir + '/hmmsearch/' + genome_id + '/' + genome_id + '.allhits.parse', sep='\t', index=False)

    return outdir + '/hmmsearch/' + genome_id + '/' + genome_id + '.allhits.parse'

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
    hmmsearchdir = os.path.join(outdir, 'hmmsearch')
    subdirectories = os.listdir(hmmsearchdir)
    subdirectories_wpath = list(map(lambda dir: os.path.join(hmmsearchdir, dir), subdirectories))

    def get_outfile(dir):
        all_files = os.listdir(dir)
        hmmoutfile = list(filter(lambda x: '_hmmsearch.out' in x, all_files))[0]
        return hmmoutfile

    hmmoutfiles = list(map(get_outfile, subdirectories_wpath))

    parsed_hmm_outfiles = list(p.map(lambda outfile: parse_hmmdomtbl(outdir, outfile, threshold, best),
                                     hmmoutfiles))
    #parsed_hmm_outfiles_multidomain = list(p.map(lambda outfile: parse_hmmdomtbl_multidomain(outdir, outfile, threshold),
    #                                 hmmoutfiles))
    return parsed_hmm_outfiles


def run_hmms(fastafile, outdir, threshold, best, cut_nc, cut_ga):

    fastaoutdir = outdir + '/hmmsearch/' + fastafile.split('/')[-1].split('.faa')[0].split('.fna')[0].split('.fasta')[0].split('.fa')[0]

    # Make outdir for HMMs
    if not os.path.exists(fastaoutdir):
        os.system('mkdir ' + fastaoutdir)
    #Make symbolic link
    if len(list(filter(lambda x: '.faa' in x, os.listdir(fastaoutdir)))) == 0:
        os.system('ln -s ' + fastafile + ' ' + fastaoutdir + '/')

    # Run all HMMs for fastafile
    return run_hmmsearch(fastafile, outdir, threshold, best, cut_nc, cut_ga)

def mark_with_threshold(all_df, hmm_thresh_dict):

    if hmm_thresh_dict is not  None:
        all_df['above_threshold'] = all_df.apply(lambda x: x.overall_bitscore >= float(hmm_thresh_dict[x.family_hmm]), axis=1)

    return all_df

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
                genome_id = genome_hit.split('|')[0]
            except:
                print(genome_id)
                print(type(genome_id))

            hits[hmm][genome_id] += 1


    #Write it to tsv in outdir without index (annoying)
    colnames = hits.columns.values.tolist()

    hits.to_csv(outdir + '/HITSTABLE.tsv', sep='\t')

def format_headers(protdir, outdir, threads):
    #Make folder for labeled protein files
    labeled_proteins = os.path.abspath(os.path.join(outdir, 'labeled_proteins'))
    os.system('mkdir ' + labeled_proteins)
    print("Protdir: ", protdir)
    print("Labeled proteins dir: ", labeled_proteins)
    def labeler(fastafile):
        new_recs = []
        genome_id = fastafile.split('/')[-1].split('.fa')[0].split('.fna')[0]
        for rec in SeqIO.parse(fastafile, 'fasta'):
                #Assume labeling has already been done; return
                if rec.id.split('|')[0] == genome_id:
                        SeqIO.write(SeqIO.parse(fastafile, 'fasta'), os.path.join(labeled_proteins, fastafile.split('/')[-1]), 'fasta')
                        return
                new_rec = rec
                new_rec.id = genome_id + '|' + rec.id
                new_recs.append(new_rec)
        SeqIO.write(new_recs, os.path.join(labeled_proteins, fastafile.split('/')[-1]), 'fasta')
        return

    p = Pool(threads)
    protfiles_wpath = [os.path.join(protdir, x) for x in os.listdir(protdir)]
    #print("TEST RUN:")
    #print("Protein file:", protfiles_wpath[0])
    #labeler(protfiles_wpath[0])
    #sys.exit()
    def is_fasta(filename):
        if filename.endswith('.fasta') or filename.endswith('.faa') or filename.endswith('.fa') or filename.endswith('.fasta') or filename.endswith('.fna'):
            return True
        else:
            return False
    label = list(p.map(labeler, filter(is_fasta, protfiles_wpath)))
    return labeled_proteins

def ensure_synteny(all_hits, outdir, threads):
    genome_ids = all_hits.genome_id.unique()
    genome_abundances = all_hits.genome_id.value_counts()
    all_hits['orfnum'] = all_hits.orf_id.apply(lambda x: int(x.split('_')[-1]))

    #Find genomes with >=8 RP hits
    abundances_past_threshold = genome_abundances[(genome_abundances >= 8.0)]
    all_hits_thresh = all_hits[all_hits.genome_id.isin(abundances_past_threshold.index.tolist())]

    #Get scaffold IDs
    all_hits_thresh['scaffold_id'] = all_hits_thresh.orf_id.apply(lambda x: '_'.join(x.split('_')[0:-1]))

    good_ids = pd.Series(all_hits_thresh.genome_id.unique().tolist())

    #How many scaffolds are there with ribosomal hits in each genome?
    num_scaffolds = []
    for good_id in good_ids.tolist():
        num_scaffolds_id = len(all_hits_thresh[all_hits_thresh.genome_id == good_id].scaffold_id.unique())
        num_scaffolds.append(num_scaffolds_id)

    #Get a dataframe with genome_id, number_of_scaffolds_with_RPs
    scaf_series = pd.Series(num_scaffolds)
    scaf_df = pd.concat([good_ids, scaf_series], axis=1)
    scaf_df.columns = ['genome_id', 'num_scaffolds']
    scaf_df['num_hits'] = scaf_df.genome_id.apply(lambda x: abundances_past_threshold[x])

    okay_genomes = scaf_df[scaf_df.num_scaffolds == 1].genome_id.tolist()

    questionable_genomes = scaf_df[scaf_df.num_scaffolds > 1]

    questionable_genome_ids = questionable_genomes.genome_id.tolist()

    scaf_counts_df = pd.DataFrame(columns=['scaffold_id', 'scaffold_count', 'genome_id'])

    for questionable_genome in questionable_genome_ids:
        #Get reduced DF for that genome
        red_df = all_hits_thresh[all_hits_thresh.genome_id == questionable_genome]
        #How many hits are on each scaffold?
        candidate_scaffold_counts = red_df.scaffold_id.value_counts()

        #Make a DF to access that data- how many RPs on each scaffold?
        cur_scaf_df = pd.concat([pd.Series(candidate_scaffold_counts.index.tolist()), pd.Series(candidate_scaffold_counts.values)], axis=1)
        cur_scaf_df.columns = ['scaffold_id', 'scaffold_count']

        #weird way of just making column of genome_id
        cur_scaf_df['genome_id'] = cur_scaf_df['scaffold_id'].apply(lambda x: questionable_genome)
        scaf_counts_df = pd.concat([scaf_counts_df, cur_scaf_df])

    reduced_scaffolds = pd.DataFrame(columns=['scaffold_id', 'scaffold_count', 'genome_id'])

    for genome in scaf_counts_df.genome_id.unique():
        red_df = scaf_counts_df[scaf_counts_df.genome_id == genome]
        red_df = red_df.sort_values(by='scaffold_count', ascending=False)
        reduced_scaffolds = reduced_scaffolds.append(red_df.iloc[0])

    reduced_scaffolds = reduced_scaffolds[reduced_scaffolds.scaffold_count >= 7]

    okay_scaffolds = scaf_df[(scaf_df.genome_id.isin(okay_genomes)) & (~scaf_df.genome_id.str.contains('yelton'))]

    all_good_hits = all_hits[all_hits.genome_id.isin(okay_genomes)]

    #Now for the more complicated part - getting the right scaffold from the other genomes

    chosen_hits = all_hits[all_hits.genome_id.isin(reduced_scaffolds.genome_id.unique())]

    chosen_hits['scaffold_id'] = chosen_hits.orf_id.apply(lambda x: '_'.join(x.split('_')[0:-1]))
    chosen_hits['orfnum'] = chosen_hits.orf_id.apply(lambda x: x.split('_')[-1]).copy()
    chosen_hits = chosen_hits[chosen_hits.scaffold_id.isin(reduced_scaffolds.scaffold_id)]



    all_good_hits = all_good_hits.append(chosen_hits.drop(['scaffold_id', 'orfnum'], axis=1))

    #Write to all_hits_evalues_df for Concatenate_And_Align.py
    all_good_hits.to_csv(os.path.join(outdir, 'all_hits_evalues_df.tsv'), sep='\t', index=False)

    return

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
    ribosomal_synteny = args.synteny

    #Make output directory
    if not os.path.exists(outdir):
        os.system('mkdir ' + outdir)
        outdir = str(Path(outdir).absolute())
    else:
        outdir = str(Path(outdir).absolute())


    have_proteins = True if prodigaldir is not None else False
    if have_proteins:
        prodigaldir = format_headers(prodigaldir, outdir, threads)

    p = Pool(threads)

    if not have_proteins:
        # Get list of paths of all fastas
        fastalist_wpath = list(map(lambda file: os.path.join(nucdir, file), os.listdir(nucdir)))
        #mutation because i'm lazy; probably will bite me in teh ass
        fastalist_wpath = list(filter(lambda x: x.endswith('.fna') or x.endswith('.fasta') or x.endswith('.fa'), fastalist_wpath))

    # Get list of all fastas
    #fastalist = list(map(lambda file: file.split('.f')[0], os.listdir(fastadir)))

    # Get list of paths of all HMM files
    hmmlist_wpath = list(map(lambda file: os.path.join(hmmdir, file), filter(lambda x: x.endswith('.hmm'), os.listdir(hmmdir))))

    # Get list of all HMMs
    hmmlist = list(map(lambda file: file.split('.hmm')[0], filter(lambda x: x.endswith('.hmm'),os.listdir(hmmdir))))

    if not already_scanned:
        #Previous workflow used hmmsearch
            #Generate binary files for hmmsearch
            #hmm_thresh_dict = hmmprfess(hmmlist_wpath, outdir, cut_nc, cut_ga)
        hmm_thresh_dict = gather_hmms(hmmlist_wpath, outdir, cut_nc, cut_ga)
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
            format_headers(os.path.join(outdir, 'proteins'), outdir, threads)


        if not have_proteins:
            protdir = outdir + '/labeled_proteins'
        else:
            protdir = os.path.abspath(prodigaldir)
            #format_headers(protdir, outdir, threads)

        protlist_wpath = list(map(lambda file: os.path.join(protdir, file), os.listdir(protdir)))
        protlist_wpath = list(filter(lambda file: file.endswith('.faa') or file.endswith('.fa') or file.endswith('.fasta'), protlist_wpath))
        #Get list of protein files without full path
        protlist = list(map(lambda path: path.split('/')[-1].split('.fna')[0].split('.fa')[0].split('.fasta')[0],
                        protlist_wpath))

        #Make directory to store hmmsearch outfiles
        if not os.path.exists(outdir + '/hmmsearch/'):
            os.system('mkdir ' + outdir + '/hmmsearch/')

        #Make sure you get rid of any Nones
        parsed_hmm_outfiles = list(filter(lambda x: x is not None,
                                    list(p.map(lambda x: run_hmms(x, outdir, threshold, best, cut_nc, cut_ga),
                                    protlist_wpath))))

        all_df_list = list(p.map(lambda x: pd.read_csv(x, sep='\t'), parsed_hmm_outfiles))
        all_df_init = pd.concat(all_df_list, sort=False)
        #print(all_df_init.orf_id.head())
        #sys.exit()
        if hmm_thresh_dict is not None:
            all_df_thresh = mark_with_threshold(all_df_init, hmm_thresh_dict)
            all_df = all_df_thresh[all_df_thresh['above_threshold']]
        else:
            all_df = all_df_init

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

    #Ensure 16RP synteny
    if ribosomal_synteny:
        ensure_synteny(all_df, outdir, threads)


    if not no_seqs:
        rec_ids_list_by_hmm = extract_hits_4(all_df, threads, protdir, outdir)
        make_hitstable_df(rec_ids_list_by_hmm, hmmlist, protlist, outdir)
    else:
        print("no_seqs specified; HITSTABLE.tsv not created.")
    print("You did it!")
    os.system('rm -rf ' + os.path.join(outdir, 'pullseq_tmp'))
    p.terminate()
    sys.exit()

if __name__ == "__main__":
    main()
