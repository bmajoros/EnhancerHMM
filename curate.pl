#!/bin/env perl

my $MUMMIE=$ENV{"MUMMIE"};

die "curate.pl <in-dir> <outfile>\n" unless @ARGV==2;
my ($indir,$outfile)=@ARGV;

open(OUT,">$outfile") || die "can't write to file: $outfile\n";
my @files=`ls $indir/*.fastb`;

my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i];
  chomp $file;
  next unless -e $file;
  my $filestem=$file;
  if($filestem=~/([^\/]+\.fastb)/) {}
  print "cat $file | $MUMMIE/fastb-to-xgraph.pl $filestem\n";
  system("cat $file | $MUMMIE/fastb-to-xgraph.pl -t $filestem");
  my $cat=<STDIN>;
  chomp $cat;
  print OUT "$cat\t$file\n";
}
close(OUT);


