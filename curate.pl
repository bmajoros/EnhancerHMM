#!/bin/env perl

my $MUMMIE=$ENV{"MUMMIE"};

die "curate.pl <likelihoods.txt> <outfile>\n" unless @ARGV==2;
my ($likelihoodFile,$outfile)=@ARGV;

open(OUT,">$outfile") || die "can't write to file: $outfile\n";
my @files;
open(IN,$likelihoodFile) || die "can't open $likelihoodFile\n";
while(<IN>) {
  chomp; my @fields=split; next unless @fields>=2;
  my ($LL,$file)=@fields;
  push @files,[$LL,$file];
}
close(IN);

@files=sort {$b->[0] <=> $a->[0]} @files;

my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i]->[1];
  next unless -e $file;
  next if $file=~/\.t00\./;
  system("cat $file | $MUMMIE/fastb-to-xgraph.pl");
  my $cat=<STDIN>;
  chomp $cat;
  print OUT "$cat\t$dir/$file\n";
}
close(OUT);


