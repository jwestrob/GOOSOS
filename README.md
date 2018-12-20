# GOOSOS
GOOSOS: Get ribOsOmal proteinS from prOtein faStas

# Why the heck did you name it that?

Ask my friend Alex T.

# How is it pronounced?

Asking the important questions as always, I see. I pronounce it like one might pronounce 'Gooses', even though that's not a real word. Just roll with it, I guess.

# What is Gooses supposed to do, then?

Given a set of protein fastas and a set of HMMs (ribosomal or not- SCGs or other stuff would be fine), this script will pull out the HMM hits from your protein fastas, and put them in a certain order so you can align them and later concatenate them to do phylogenetic analysis.

This script also yields a 'hits table', which is a matrix (in .tsv format) of hits per genome. This is so you can quality-check later and remove genomes which don't have enough hits to be included in your phylogenetic analysis.

I hope this will be useful for other things as well, though, and if you see things that you'd like this script to do which aren't already incorporated, just let me know and I'll be happy to help out if it's not gonna take, like, a week.

# Tell me about the other cool fun things that GOOSOS can do!

So glad to see the enthusiasm. Here's a list. I know you love those.

- Evalue thresholds (-evalue): This is the default evalue cutoff for HMMsearch. I recommend something like 1e-10 or lower. (HMMsearch accepts scientific notation, so type it out like that. If you hate conformity I guess you can type out 0.0000000001)

- Extract Best Hit (-best_only) [DEFAULT] : Pulls out the most relevant hit (i.e., the hit with the lowest e-value) per HMM per genome. That is, if a genome has two hits for a particular HMM, it takes the hit with the lower evalue.

This list will be updated iteratively as I put more stuff in here.

# Wow! Cool! I want to use this in my paper!

Sick yo, but please cite this repository. I'll be putting the GPL license on it in a little bit once I actually get the thing running.
