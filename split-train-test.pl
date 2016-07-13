#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <indir> <num-test> <outdir-train> <outdir-test>\n" unless @ARGV==4;
my ($indir,$numTest,$outTrain,$outTest)=@ARGV;

my @files=`ls $indir`;
my $n=@files;
for(my $i=0 ; $i<$n ; ++$i) {
  my $file=$files[$i]; chomp $file;
  if($i<$numTest) { system("ln -s $indir/$file $outTest/$file") }
  else { system("ln -s $indir/$file $outTrain/$file") }
}





