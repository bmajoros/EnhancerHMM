#!/usr/bin/perl
use strict;
use ProgramName;

my $MUMMIE=$ENV{"MUMMIE"};

my $name=ProgramName::get();
die "$name <from-dir> <to-dir> <pos.hmm> <neg.hmm>\n" unless @ARGV==4;
my ($fromDir,$toDir,$posHMM,$negHMM)=@ARGV;

process($fromDir,$toDir,$posHMM,$negHMM,1);

#########################################################################
sub process
{
  my ($fromDir,$toDir,$posHMM,$negHMM,$label)=@_;
  my @files=`ls $fromDir`;
  my $n=@files;
  for(my $i=0 ; $i<$n ; ++$i) {
    my $file=$files[$i]; chomp $file;
    next unless $file=~/\.fastb$/;
    my $numer=`$MUMMIE/get-likelihood $posHMM $dir/$file`;
    my $denom=`$MUMMIE/get-likelihood $negHMM $dir/$file`;
    my $ratio=$numer-$denom;
    print "$ratio\n";
  }
}








