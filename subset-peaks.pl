#!/bin/env perl
use strict;

die "subset-peaks.pl <file-list.txt> <num-subsets> <outdir>\n" unless @ARGV==3;
my ($filelist,$numSets,$outDir)=@ARGV;

my @files;
open(IN,$filelist) || die $filelist;
while(<IN>) {
  chomp; next unless(/\S+/);
  push @files,$_;
}
close(IN);

my $numFiles=@files;
my $filesPerSet=$numFiles/$numSets;
my $nextSet=1;
my $file=">$outDir/files$nextSet.txt";
open(OUT,$file) || die $file;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  if($i>$nextSet*$filesPerSet) {
    close(OUT);
    ++$nextSet;
    my $file=">$outDir/files$nextSet.txt";
    open(OUT,$file) || die $file;
  }
  my $file=$files[$i];
  print OUT "$file\n";
}
close(OUT);


