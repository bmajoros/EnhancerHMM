#!/usr/bin/perl
use strict;
use ProgramName;
use TempFilename;

my $MARKOV="/data/reddylab/software/meme_4.10.0/bin/fasta-get-markov";
my $MAST="/data/reddylab/software/meme_4.10.0/bin/mast";
my $EXTRACT="/home/bmajoros/MUMMIE/fastb-to-fasta.pl";
my $TEMP_FASTA=TempFilename::generate();
my $TEMP_BACKGROUND=TempFilename::generate();
my $TEMP_OUT=TempFilename::generate();

my $name=ProgramName::get();
die "$name <in-fastb-dir> <out-fastb-dir>" unless @ARGV==2;
my ($inDir,$outDir)=@ARGV;

my @files=`ls $inDir`;
my $n=@files;
for(my $i=0 ; $i<$n ; ++$i) {
  my $file=$files[$i]; chomp $file;
  next unless $file=~/\.fastb$/;
  system("$EXTRACT $inDir/$file $TEMP_FASTA");

  system("$MARKOV -m 1 $TEMP_FASTA > $TEMP_BACKGROUND");

  my $MOTIF="main-motifs/GR.meme";
  system("$MAST $MOTIF $TEMP_FASTA -hit_list -mt 0.05 -bfile $TEMP_BACKGROUND > $TEMP_OUT");

  exit;
}

unlink($TEMP_FASTA); unlink($TEMP_BACKGROUND); unlink($TEMP_OUT);


