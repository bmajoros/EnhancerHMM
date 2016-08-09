#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in-dir> <out-dir>\n" unless @ARGV==2;
my ($indir,$outdir)=@ARGV;

my %times;
my @files=`ls $indir`;
my $n=@files;
for(my $i=0 ; $i<$n ; ++$i) {
  my $file=$files[$i]; chomp $file;
  next unless $file=~/\.fastb$/;
  $file=~/(\S+)\.t(\d+).standardized.fastb/ || die $file;
  my ($id,$time)=($1,$2);
  my $prev=$times{$id};
  if(!defined($prev) || $prev<$time) { $times{$id}=$time }
}

my @IDs=keys %times;
my $n=@IDs;
for(my $i=0 ; $i<$n ; ++$i) {
  my $id=$IDs[$i];
  my $time=$times{$id};
  my $file="$id.t$time.standardized.fastb";
  system("cp $indir/$file $outdir\n");
}




