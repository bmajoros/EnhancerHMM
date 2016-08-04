#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in-fastb-dir> <out-fastb-dir>" unless @ARGV==2;
my ($inDir,$outDir)=@ARGV;





"/data/reddylab/software/meme_4.10.0/bin/fasta-get-markov -m 1 fastafile.fa > background.out"

"/data/reddylab/software/meme_4.10.0/bin/mast GR.motif fastafile.fa -hit_list -mt 0.05 -bfile background.out  > mast_results.out"



