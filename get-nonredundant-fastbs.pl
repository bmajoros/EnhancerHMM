#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in-dir> <out-dir>\n" unless @ARGV==2;
my ($indir,$outdir)=@ARGV;

my %seen;
my @files=`ls $indir`;
my $n=@files;
for(my $i=0 ; $i<$n ; ++$i) {
  my $file=$files[$i]; chomp $file;
  next unless $file=~/\.fastb$/;
  $file=~/(\S+)\.t\d+.standardized.fastb/ || die $file;
  my $id=$1;
  next if $seen{$id};
  system("cp $indir/$file $outdir");
  $seen{$id}=1;
}





