#!/bin/env perl
use strict;
use ProgramName;

my $BASE="/home/bmajoros/GGR/p300";

my $name=ProgramName::get();
die "$name <trajectories.trxt>\n" unless @ARGV==1;
my ($infile)=@ARGV;

open(IN,$infile) || die "can't open $infile\n";
while(<IN>) {
  chomp; my @fields=split; next unless @fields>=13;
  my $peak=shift @fields;
  my $class=classify(\@fields);
  print "$peak\t$class\n";
}
close(IN);



# =====================================================================
sub classify {
  my ($array)=@_;
  my $firstHalf=$count($array,0,6);
  my $secondHalf=$count($array,6,12);
  my $total=$count($array,0,12);
  if($total>=11) { return "on" }
  if($total<=1) { return "off" }
  if($firstHalf<=2 && $secondHalf>=4) { return "increasing" }
  if($firstHalf>=4 && $secondHalf<=2) { return "decreasing" }
  return "ambiguous";
}
# =====================================================================
sub count {
  my ($array,$begin,$end)=@_;
  my $count=0;
  for(my $i=$begin ; $i<$end ; ++$i) {
    if($array->[$i]>0.5) { ++$count }
  }
  return $count;
}
# =====================================================================
# =====================================================================
# =====================================================================
# =====================================================================
# =====================================================================
# =====================================================================

