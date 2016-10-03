#!/usr/bin/perl
use strict;
$|=1;

my $BASE="/home/bmajoros/GGR/p300";
my $TRAJECTORIES="$BASE/trajectories-1.txt";

my @times=(0,0.5,1,2,3,4,5,6,7,8,10,12);

my (@counts,$sampleSize);
open(IN,$TRAJECTORIES) || die $TRAJECTORIES;
while(<IN>) {
  chomp; my @fields=split; next unless @fields>=13;
  shift @fields;
  my ($hasZero,$hasOne);
  for(my $i=0 ; $i<12 ; ++$i) {
    my $x=$fields[$i];
    if($x==1) { $hasOne=1 } elsif($x==0) { $hasZero=1 } }
  next unless $hasZero && $hasOne;
  for(my $i=0 ; $i<12 ; ++$i) { $counts[$i]+=$fields[$i] }
  ++$sampleSize;
}
close(IN);

for(my $i=0 ; $i<12 ; ++$i) {
  my $p=$counts[$i]/$sampleSize;
  $p=int($p*100+5/9)/100;
  my $time=$times[$i];
  print "$time\t$p\n";
}
print STDERR "sample size $sampleSize\n";




