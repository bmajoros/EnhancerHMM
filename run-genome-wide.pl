#!/usr/bin/perl
use strict;
use ProgramName;
$|=1;

my $MUMMIE=$ENV{"MUMMIE"};

my $name=ProgramName::get();
die "$name <dir> <pos.hmm> <neg.hmm>\n" unless @ARGV==3;
my ($dir,$posHMM,$negHMM)=@ARGV;

process($dir,$posHMM,$negHMM);

#########################################################################
sub process
{
  my ($dir,$posHMM,$negHMM)=@_;
  my @files=`ls $dir`;
  my $n=@files;
  for(my $i=0 ; $i<$n ; ++$i) {
    my $file=$files[$i]; chomp $file;
    next unless $file=~/(\S+)\.fastb$/;
    my $id=$1;
    my $numer=`$MUMMIE/get-likelihood $posHMM $dir/$file`;
    my $denom=`$MUMMIE/get-likelihood $negHMM $dir/$file`;
    my $ratio=$numer-$denom;
    print "$id\t$ratio\n";
  }
}








