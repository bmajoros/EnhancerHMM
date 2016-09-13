#!/usr/bin/perl
use strict;
use Fastb;

my $BASE="/home/bmajoros/GGR/p300";
my $POS="$BASE/motifs-pos";
my $NEG="$BASE/motifs-neg";
my $RAW="$BASE/raw";
my $POS_OUT="$BASE/pos-0hr-3hr";
my $NEG_OUT="$BASE/neg-0hr-3hr";

process($POS,$POS_OUT);
process($NEG,$NEG_OUT);

sub process {
  my ($inDir,$outDir)=@_;
  my @files=`ls $inDir`;
  my $numFiles=@files;
  for(my $i=0 ; $i<$numFiles ; ++$i) {
    my $file=$files[$i]; chomp $file;
    next unless $file=~/(iter\d+_peak\d+\.t)(\d+)\.fastb/;
    my ($stem,$time)=($1,$2);
    my $zeroFile="$RAW/${stem}00.fastb";
    my $threeFile="$inDir/${stem}3.fastb";
    my $zeroFastb=new Fastb($zeroFile);
    my $threeFastb=new Fastb($threeFile);
    $threeFastb->dropTrack("dna");
    my $dnase0=$zeroFastb->getTrackByName("DNase");
    my $p300_0=$zeroFastb->getTrackByName("EP300");
    $dnase0->rename("dnase0");
    $p300_0->rename("p300_0");
    $threeFastb->addTrack($dnase0);
    $threeFastb->addTrack($p300_0);
    my $outfile="$outDir/$file";
    $threeFastb->save($outfile);
  }
}


