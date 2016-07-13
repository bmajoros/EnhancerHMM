#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <indir> <num-test> <outdir-train> <outdir-test>\n" unless @ARGV==4;
my ($indir,$numTest,$outTrain,$outTest)=@ARGV;

if($indir=~/^\//) { die "please use a relative path, not absolute\n" }

my @files=`ls $indir`;
my $n=@files;
for(my $i=0 ; $i<$n ; ++$i) {
  my $file=$files[$i]; chomp $file;
  if($i<$numTest) { system("cd $outTest ; ln -s ../$indir/$file") }
  else { system("cd $outTrain ; ln -s ../$indir/$file") }
}





