#!/bin/env perl
use strict;
use ProgramName;
$|=1;

my $BASE="/home/bmajoros/GGR/p300";

my $name=ProgramName::get();
die "$name <trajectories.txt>\n" unless @ARGV==1;
my ($infile)=@ARGV;

open(IN,$infile) || die "can't open $infile\n";
while(<IN>) {
  chomp; my @fields=split; next unless @fields>=13;
  my $peak=shift @fields;
  my $class=classify(\@fields);
  #print "$peak\t$class\n";
  print "$peak\t$class\t@fields\n";
}
close(IN);



# =====================================================================
sub classify {
  my ($array)=@_;
  my $string=join("",@$array);
  if($string=~/^(0+)(1+)/) {
    my $leftLen=length($1); my $rightLen=length($2);
    if($leftLen>=2 && $rightLen>=4) { return "increasing" }
    else { return "ambiguous" } }
  if($string=~/^(1+)(0+)/) {
    my $leftLen=length($1); my $rightLen=length($2);
    if($leftLen>=2 && $rightLen>=4) { return "decreasing" }
    else { return "ambiguous" } }
  elsif($string=~/^0+$/) { return "off" }
  elsif($string=~/^1+$/) { return "on" }
  else { return "ambiguous" }
}
# =====================================================================
sub classify_old {
  my ($array)=@_;
  my $firstHalf=count($array,0,6);
  my $secondHalf=count($array,6,12);
  my $total=count($array,0,12);
  if($total>=11) { return "on" }
  if($total<=1) { return "off" }
#  if($firstHalf<=1 && $secondHalf>=5) { return "increasing" }
#  if($firstHalf>=5 && $secondHalf<=1) { return "decreasing" }
  if($firstHalf-$secondHalf>=3) { return "decreasing" }
  elsif($secondHalf-$firstHalf>=3) { return "increasing" }
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

