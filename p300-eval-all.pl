#!/usr/bin/perl
use strict;
$|=1;

my $BASE="/home/bmajoros/GGR/p300";
my $LS="$BASE/ls-raw.txt";
my $HMM="$BASE/model/trained-pos5.hmm";
my $MUMMIE=$ENV{"MUMMIE"};

my @files=`ls no-dna`;
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i];
  chomp $file; next unless($file=~/\.fastb$/);
  my $file="no-dna/$file";
#die "$MUMMIE/get-likelihood $HMM $file";
  my $LL=0+`$MUMMIE/get-likelihood $HMM $file`;
  print "$LL\t$file\n";
}


