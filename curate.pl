#!/bin/env perl

my $MUMMIE=$ENV{"MUMMIE"};

die "curate.pl <fastb-dir> <outfile>\n" unless @ARGV==2;
my ($dir,$outfile)=@ARGV;

open(OUT,">$outfile") || die "can't write to file: $outfile\n";
my @files=`ls $dir`;
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i]; chomp $file;
  next unless $file=~/\.fastb$/;
  system("cat $dir/$file | $MUMMIE/fastb-to-xgraph.pl");
  my $cat=<STDIN>;
  chomp $cat;
  print OUT "$cat\t$dir/$file\n";
}
close(OUT);

