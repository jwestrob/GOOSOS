# GOOSOS
GOOSOS: Get ribOsOmal proteinS from prOtein faStas

# Why the heck did you name it that?

Ask my friend Alex T.

# How is it pronounced?

Asking the important questions as always, I see. I pronounce it like one might pronounce 'Gooses', even though that's not a real word. Just roll with it, I guess.

# What is Gooses supposed to do, then?

Given a set of nucleotide fastas and a set of HMMs (ribosomal or not- SCGs or other stuff would be fine), this script will predict protein sequences for your nucleotide sequences, pull out the HMM hits from the resulting protein fastas, and put them in a certain order so you can align them and later concatenate them to do phylogenetic analysis.

GOOSOS compresses the HMMs you provide into a binary database, then searches the predicted protein sequences for hits to those HMMs using HMMscan. I then parse out all those

This script yields two tab-separated dataframe-style output files: `all_hits_evalues_df.tsv` and `HITSTABLE.tsv`.

- `all_hits_evalues_df.tsv` has a list of hits organized by genome of origin, along with ORF numbers indicating the contig/scaffold they originated from in the original nucleotide data, so that you can check for synteny among your hits. This is intended for ribosomal protein analyses, but might hopefully be useful when looking at other proteins too. (Trying to find an operon, perhaps?) This dataframe also contains helpful e-value information on your hits, so that you can determine a sensible cutoff value if you haven't yet determined one (or don't want to just pick a threshold before looking at your data).

    - Wondering what all this nonsense about domain 1 and 2 is? Well, lots of hits (at least for the ribosomal sequences I've been looking at) have two different domains, so I've included information on both of those domains when relevant. I am not particularly interested in cases where you might have more than two domain hits for a single ORF, but if you are and you're concerned about this let me know I guess. I can adapt this dataframe to include potential relevant information. If you want to check whether there are ORFs with >2 hits, look at the HMMscan outfiles (output/hmmscan/).

- `HITSTABLE.tsv` indicates how many hits for each HMM are in each genome. I use this for quality control once I've done all my scanning- for instance, you might not want to include genomes that have 2 distinct hits for each of your ribosomal HMMs when doing phylogenomic analyses. (Looking at you, SAR202 megabins)

I hope this will be useful for other things as well, though, and if you see things that you'd like this script to do which aren't already incorporated, just let me know and I'll be happy to help out if it's not gonna take, like, a week.

# Tell me about the other cool fun things that GOOSOS can do!

So glad to see the enthusiasm. Here's a list. I know you love those.

- Evalue thresholds (-evalue): This is the default evalue cutoff for HMMsearch. I recommend something like 1e-10 or lower. (HMMsearch accepts scientific notation, so type it out like that. If you hate conformity I guess you can type out 0.0000000001)

- Extract Best Hit (-best) [DEFAULT] : Pulls out the most relevant hit (i.e., the hit with the lowest e-value) per HMM per genome. That is, if a genome has two hits for a particular HMM, it takes the hit with the lower evalue.

- If you already scanned with HMMsearch (-already_scanned) : Don't scan again! Skip to the next part, using the same output directory.

- If you don't want to write seqs and only want a hits table (-no_seqs) : Saves time if you're trying to calibrate your e-value threshold.

This list will be updated iteratively as I put more stuff in here.

# But there's another program in here! What the heck dude?

This additional program (`Concatenate_And_Align.py`) takes the fastas you generated with GOOSOS, filters out genomes that don't have enough sequences (as determined by the `-threshold` flag), aligns them, sorts them, concatenates the alignments, and yields a partition file (`partitions.nex`) that you can use when building your phylogenetic tree.



# What are the dependencies?

Well, python obviously. I recommend using pip to install all this (which can be itself used through Anaconda/Miniconda).

The packages you will need to install (using a variant of the command 'pip install [PACKAGE]') are:

- pathos (Enables multithreading)
- Biopython
- Pandas
- Argparse
- Mafft (If you want to align sequences)

Everything else should come standard with your install of python.

# Wow! Cool! I want to use this in my paper!

Sick yo, but please cite this repository. I'll be putting the GPL license on it in a little bit once I actually get the thing running.
