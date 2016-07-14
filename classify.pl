#!/usr/bin/perl
use strict;
use ProgramName;

my $ROC="tmp.roc";
my $MUMMIE=$ENV{"MUMMIE"};

my $name=ProgramName::get();
die "$name <test-pos-dir> <test-neg-dir> <pos.hmm> <neg.hmm>\n" unless @ARGV==4;
my ($posDir,$negDir,$posHMM,$negHMM)=@ARGV;

open(ROC,">$ROC") || die "can't write to file: $ROC\n";
process($posDir,$posHMM,$negHMM,1);
process($negDir,$posHMM,$negHMM,0);
close(ROC);

#########################################################################
sub process
{
  my ($dir,$posHMM,$negHMM,$label)=@_;
  my @files=`ls $dir`;
  my $n=@files;
  for(my $i=0 ; $i<$n ; ++$i) {
    my $file=$files[$i]; chomp $file;
    next unless $file=~/\.fastb$/;
    my $numer=`$MUMMIE/get-likelihood $posHMM $dir/$file`;
    my $denom=`$MUMMIE/get-likelihood $negHMM $dir/$file`;
    my $ratio=$numer-$denom;
    print ROC "$ratio\t$label\n";
  }
}








