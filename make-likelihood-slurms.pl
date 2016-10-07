#!/usr/bin/perl
use strict;
use SlurmWriter;

my $BASE="/home/bmajoros/GGR/p300";
my $SLURM_DIR="$BASE/slurms/likelihood-slurms";
my $PROGRAM="$BASE/src/get-likelihoods.pl";
my $FASTB="$BASE/no-dna";
my $HMM="$BASE/model/bootstrap-neg1.hmm";

my $writer=new SlurmWriter();
my @lists=`ls $BASE/subsets`;
foreach my $list (@lists) {
  chomp $list;
  next unless $list=~/^files(\d+)\.txt$/;
  my $ID=$1;
  my $dir="$BASE/likelihoods";
  my $infile="$BASE/subsets/$list";
  my $outfile="$BASE/likelihoods/$ID.txt";
  $writer->addCommand("$PROGRAM $HMM $infile $FASTB > $outfile");
}
#$writer->mem(5000);
$writer->setQueue("new,all");
$writer->writeArrayScript($SLURM_DIR,"LIKELIHOOD",$SLURM_DIR,500);



