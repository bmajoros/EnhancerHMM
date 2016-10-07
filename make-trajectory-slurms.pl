#!/usr/bin/perl
use strict;
use SlurmWriter;

my $BASE="/home/bmajoros/GGR/p300";
my $SLURM_DIR="$BASE/slurms/trajectory-slurms";
my $PROGRAM="$BASE/src/p300-trajectory-posteriors.pl";
my $FASTB="$BASE/no-dna";
my $POS_HMM="$BASE/model/bootstrap-pos5.hmm";
#my $NEG_HMM="$BASE/model/bootstrap-neg5.hmm";
my $NEG_HMM="$BASE/model/bootstrap-neg1.hmm";

my $writer=new SlurmWriter();
my @lists=`ls $BASE/subsets`;
foreach my $list (@lists) {
  chomp $list;
  next unless $list=~/^files(\d+)\.txt$/;
  my $ID=$1;
  my $infile="$BASE/peak-subsets/$list";
  my $outfile="$BASE/trajectories/$ID.txt";
  $writer->addCommand("$PROGRAM $POS_HMM $NEG_HMM $infile $FASTB > $outfile");
}
#$writer->mem(5000);
$writer->setQueue("new,all");
$writer->writeArrayScript($SLURM_DIR,"TRAJECTORIES",$SLURM_DIR,500);



