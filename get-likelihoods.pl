#!/bin/env perl
use strict;

my $MUMMIE=$ENV{"MUMMIE"};
my $PROGRAM="$MUMMIE/get-likelihood";

die "get-likelihoods.pl <*.hmm> <list.txt> <dir>\n" unless @ARGV==3;
my ($HMM,$listFile,$DIR)=@ARGV;

open(IN,$listFile) || die "can't open $listFile\n";
while(<IN>) {
  chomp; next unless(/\.fastb$/);
  my $file=$_;
  my $path="$DIR/$file";
  my $LL=0+`$PROGRAM $HMM $path`;
  print "$LL\t$file\n";
}
close(IN);



